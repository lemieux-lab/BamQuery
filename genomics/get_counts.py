import warnings
warnings.filterwarnings("ignore")
import os, logging, time, subprocess, pickle, multiprocessing, os, _thread, csv, collections, pysam, copy
import genomics.get_alignments as get_alig
import pandas as pd
from pathos.multiprocessing import ProcessPool
import utils.useful_functions as uf

__author__ = "Maria Virginia Ruiz Cuevas"
__email__ = "maria.virginia.ruiz.cuevas@umontreal.ca"


NUM_WORKERS =  multiprocessing.cpu_count()

path_to_lib = '/'.join(os.path.abspath(__file__).split('/')[:-3])+'/lib/'

class GetCounts:

	def __init__(self, path_to_output_folder, name_exp, mode, light, peptides_by_type):
		self.light = light

		if self.light:
			self.path_to_output_temps_folder = path_to_output_folder+'res_light/temps_files/'
		else:
			self.path_to_output_temps_folder = path_to_output_folder+'res/temps_files/'
		
		self.name_exp = name_exp
		self.mode = mode
		self.path_to_output_folder_alignments = path_to_output_folder+'alignments/'
		self.peptides_by_type = peptides_by_type
		

	def ribo_counts(self, perfect_alignments, bam_files_list):

		df_counts = pd.DataFrame()
		exists = os.path.exists(self.path_to_output_temps_folder+self.name_exp+'_ribo_count.csv')
		
		order = []

		if not exists:
			t_0 = time.time()
			logging.info('Alignments on ribosome profiling information ')

			to_write = {} 
			
			perfect_alignments_to_return = {} 

			data = {}
			info_bams = []
			bams = []

			list_bam_files_order = []
			for name_sample, info_bam in sorted(bam_files_list.items(), key=lambda e: e[1][-1], reverse=False):
				info_bams.append((name_sample,info_bam))
				list_bam_files_order.append(name_sample)
			
			def get_index(sample):
				index = list_bam_files_order.index(sample)
				return index + 1

			total_samples = len(list_bam_files_order)

			keys = perfect_alignments.keys()
			
			logging.info('Total MCS mapped : %s ', str(len(keys)))

			modif_dic = {}
			for key in keys:
				split_key = key.split('_')
				peptide = split_key[0]
				position = split_key[1]
				seq = split_key[2]

				strand = perfect_alignments[key][0]
				new_key = peptide+'_'+position+'_'+strand
				
				try:
					modif_dic[new_key].append(seq)
				except KeyError:
					modif_dic[new_key] = [seq]

			keys = modif_dic.keys()
			values = modif_dic.values()

			logging.info('Total unique regions : %s ', str(len(keys)))
			
			pool = ProcessPool(nodes=NUM_WORKERS)

			for idx, bam_file in enumerate(info_bams):
				t0_bam_file = time.time()
				bams = [bam_file] * len(keys)
				results = pool.map(self.get_counts_sample, bams, keys, values)
				
				for res in results:

					for index, count_align in enumerate(res):

						if index == 0:
							peptide = count_align[0]
							alignment = count_align[1]
							sample = count_align[2]
							strand = count_align[3]
							to_print = [peptide, alignment, sample]
						else:
							count = count_align[0]
							sequence = count_align[1]
							to_print_aux = copy.deepcopy(to_print)
							to_print_aux.append(count)
							index = get_index(sample)

							key_aux = alignment+'_'+sequence

							try:
								peptide_info_to_write = to_write[peptide]

								try:
									peptide_info_to_write[key_aux][index] = count

								except KeyError:
									counts = [0]*total_samples
									to_add = [strand]
									to_add.extend(counts)
									to_add[index] = count
									peptide_info_to_write[key_aux] = to_add

							except KeyError:
								counts = [0]*total_samples
								to_add = [strand]
								to_add.extend(counts)
								to_add[index] = count
								to_write[peptide] = {key_aux: to_add}

							try:
								data[peptide].append(to_print_aux)
							except KeyError:
								data[peptide] = [to_print_aux]

							key = peptide+'_'+alignment+'_'+sequence
							if len(perfect_alignments[key][-1]) == 0:
								perfect_alignments[key][-1] = [0]*total_samples

							perfect_alignments[key][-1][index-1] = count
							perfect_alignments_to_return[key] = perfect_alignments[key]

				t1_bam_file = time.time()
				time_final = (t1_bam_file-t0_bam_file)/60.0
				logging.info('Processed Bam File : %s %s . Time : %f min', str(idx), bam_file[0], time_final)

			pool.close()
			pool.join()
			pool.clear()
			
			header = ['Peptide', 'Position', 'Strand']
			header.extend(list_bam_files_order)
			to_write_list = [header]

			for peptide, info_peptide in to_write.items():
				to_add = []
				for alignment, info_counts in info_peptide.items():
					position = alignment.split('_')[0]
					MCS = alignment.split('_')[1]
					to_add =  [peptide, position, MCS]
					to_add.extend(info_counts)
					to_write_list.append(to_add)

			df_alignments = pd.DataFrame(to_write_list[1:], columns=to_write_list[0])

			if self.light:
				name_path_info = self.path_to_output_folder_alignments+'Alignments_information_light.dic'
				name_path_info_ribo = self.path_to_output_folder_alignments+'Alignments_ribo_information_light.dic'
			else:
				name_path_info = self.path_to_output_folder_alignments+'Alignments_information.dic'
				name_path_info_ribo = self.path_to_output_folder_alignments+'Alignments_ribo_information.dic'

			with open(name_path, 'wb') as handle:
				pickle.dump(perfect_alignments, handle, protocol=pickle.HIGHEST_PROTOCOL)

			with open(name_path, 'wb') as handle:
				pickle.dump(perfect_alignments_to_return, handle, protocol=pickle.HIGHEST_PROTOCOL)

			data_ordered = []
			order_pep = []
			peptides_order = []

			for peptide_type, peptides in self.peptides_by_type.items():
				for peptide in peptides:
					try:
						data_ordered.extend(data[peptide])
						peptides_order.append(peptide)
						if peptide_type not in order_pep:
							order_pep.append(peptide_type)
							order.append(1)
						else:
							ind = order_pep.index(peptide_type)
							order[ind] += 1
					except KeyError:
						pass

			if len(data_ordered) > 0 :
				df_counts = pd.DataFrame(data_ordered, columns=['Peptides', 'Alignments', 'BAM Files', 'Read Counts'])
				df_counts = df_counts.groupby(['Peptides', 'BAM Files'])['Read Counts'].sum().reset_index()
				df_counts = df_counts.pivot("Peptides", "BAM Files", "Read Counts").reset_index()

				sorterIndex = dict(zip(peptides_order, range(len(peptides_order))))
				df_counts['order'] = df_counts_filtered['Peptides'].map(sorterIndex)
				df_counts.sort_values(['order'], ascending = [True], inplace = True)
				df_counts.drop('order', 1, inplace = True)
				df_counts.set_index('Peptides',inplace=True)
				_thread.start_new_thread(self.save_info_counts, (df_counts, to_write_list, '_ribo_count.csv'))
			
			t_2 = time.time()
			total = t_2-t_0
			logging.info('Total time run function get_counts to end : %f min', (total/60.0))
		else:
			logging.info('Count information already collected in the output folder : %s --> Skipping this step!', self.path_to_output_temps_folder+self.name_exp+'_ribo_count.csv')
			df_counts = pd.read_csv(self.path_to_output_temps_folder+self.name_exp+'_ribo_count.csv', index_col=0)

			if self.light:
				name_path = self.path_to_output_folder_alignments+'Alignments_information_light.dic'
			else:
				name_path = self.path_to_output_folder_alignments+'Alignments_information.dic'

			with open(name_path, 'rb') as fp:
				perfect_alignments = pickle.load(fp)
			
			df_alignments = pd.read_csv(self.path_to_output_temps_folder+self.name_exp+'_ribo_count_All_alignments.csv', index_col=0)

		return perfect_alignments, df_counts, order, df_alignments


	def get_counts(self, perfect_alignments, bam_files_list):

		df_counts = pd.DataFrame()
		df_counts_filtered = pd.DataFrame()

		exist_rna = os.path.exists(self.path_to_output_temps_folder+self.name_exp+'_rna_count.csv')
		exist_ribo = True

		if self.mode == 'translation':

			if self.light:
				name_path = self.path_to_output_folder_alignments+'Alignments_ribo_information.dic'
			else:
				name_path = self.path_to_output_folder_alignments+'Alignments_ribo_information_light.dic'

			with open(name_path, 'rb') as fp:
				perfect_alignments_ribo = pickle.load(fp)

			data_filtered = {}
			exist_ribo = os.path.exists(self.path_to_output_temps_folder+self.name_exp+'_rna_ribo_count.csv')

		order_f = []
		order = []

		if not exist_rna and exist_ribo:
			t_0 = time.time()
			to_write = {} 
			
			data = {}
			
			info_bams = []
			bams = []
			
			list_bam_files_order = []
			for name_sample, info_bam in sorted(bam_files_list.items(), key=lambda e: e[1][-1], reverse=False):
				info_bams.append((name_sample,info_bam))
				list_bam_files_order.append(name_sample)
			
			def get_index(sample):
				index = list_bam_files_order.index(sample)
				return index + 1

			total_samples = len(list_bam_files_order)

			keys = perfect_alignments.keys()
			
			logging.info('Total MCS mapped : %s ', str(len(keys)))

			modif_dic = {}
			for key in keys:
				split_key = key.split('_')
				peptide = split_key[0]
				position = split_key[1]
				seq = split_key[2]

				strand = perfect_alignments[key][0]
				new_key = peptide+'_'+position+'_'+strand
				
				try:
					modif_dic[new_key].append(seq)
				except KeyError:
					modif_dic[new_key] = [seq]

			keys = modif_dic.keys()
			values = modif_dic.values()

			logging.info('Total unique regions : %s ', str(len(keys)))

			pool = ProcessPool(nodes=NUM_WORKERS)

			for idx, bam_file in enumerate(info_bams):
				t0_bam_file = time.time()
				bams = [bam_file] * len(keys)
				results = pool.map(self.get_counts_sample, bams, keys, values)
				
				for res in results:

					for index, count_align in enumerate(res):

						if index == 0: 
							peptide = count_align[0]
							alignment = count_align[1]
							sample = count_align[2]
							strand = count_align[3]
							to_print = [peptide, alignment, sample]
						else:
							count = count_align[0]
							sequence = count_align[1]
							to_print_aux = copy.deepcopy(to_print)
							to_print_aux.append(count)
							
							index = get_index(sample) 

							key_aux = alignment+'_'+sequence
							try:
								peptide_info_to_write = to_write[peptide]

								try:
									peptide_info_to_write[key_aux][index] = count
								except KeyError:
									counts = [0]*total_samples
									to_add = [strand]
									to_add.extend(counts)
									to_add[index] = count
									peptide_info_to_write[key_aux] = to_add

							except KeyError:
								counts = [0]*total_samples
								to_add = [strand]
								to_add.extend(counts)
								to_add[index]  = count
								to_write[peptide] = {key_aux: to_add}

							try:
								data[peptide].append(to_print_aux)
							except KeyError:
								data[peptide] = [to_print_aux]

							key = peptide+'_'+alignment+'_'+sequence
							if len(perfect_alignments[key][-2]) == 0:
								perfect_alignments[key][-2] = [0]*total_samples

							perfect_alignments[key][-2][index-1] = count

							if self.mode == 'translation':
								try:
									count_ribo = perfect_alignments_ribo[key][-1][index-1]
									if count_ribo > 0:
										try:
											data_filtered.append(to_print_aux)
										except KeyError:
											data_filtered[peptide] = [to_print_aux]

								except KeyError:
									pass

				t1_bam_file = time.time()
				time_final = (t1_bam_file-t0_bam_file)/60.0
				logging.info('Processed Bam File : %s %s . Time : %f min', str(idx), bam_file[0], time_final)
			
			pool.close()
			pool.join()
			pool.clear()

			header = ['Peptide', 'Position', 'MCS', 'Strand']
			header.extend(list_bam_files_order)
			to_write_list = [header]
			
			for peptide, info_peptide in to_write.items():
				to_add = []
				for alignment, info_counts in info_peptide.items():
					position = alignment.split('_')[0]
					MCS = alignment.split('_')[1]
					to_add =  [peptide, position, MCS]
					to_add.extend(info_counts)
					to_write_list.append(to_add)

			df_alignments = pd.DataFrame(to_write_list[1:], columns=to_write_list[0])

			if self.light:
				name_path = self.path_to_output_folder_alignments+'Alignments_information_light.dic'
			else:
				name_path = self.path_to_output_folder_alignments+'Alignments_information.dic'

			with open(name_path, 'wb') as handle:
				pickle.dump(perfect_alignments, handle, protocol=pickle.HIGHEST_PROTOCOL)

			data_ordered = []
			order_pep = []
			peptides_order = []
			
			for peptide_type, peptides in self.peptides_by_type.items():
				for peptide in peptides:
					try:
						data_ordered.extend(data[peptide])
						peptides_order.append(peptide)
						if peptide_type not in order_pep:
							order_pep.append(peptide_type)
							order.append(1)
						else:
							ind = order_pep.index(peptide_type)
							order[ind] += 1
					except KeyError:
						pass

			if len(data_ordered) > 0 :
				# https://stackoverflow.com/questions/23482668/sorting-by-a-custom-list-in-pandas

				df_counts = pd.DataFrame(data_ordered, columns=['Peptides', 'Alignments', 'BAM Files', 'Read Counts'])
				df_counts = df_counts.groupby(['Peptides', 'BAM Files'])['Read Counts'].sum().reset_index()
				df_counts = df_counts.pivot("Peptides", "BAM Files", "Read Counts").reset_index() 

				sorterIndex = dict(zip(peptides_order, range(len(peptides_order))))
				df_counts['order'] = df_counts['Peptides'].map(sorterIndex)
				df_counts.sort_values(['order'], ascending = [True], inplace = True)
				df_counts.drop('order', 1, inplace = True)
				df_counts.set_index('Peptides',inplace=True)

				_thread.start_new_thread(self.save_info_counts, (df_counts, to_write_list, '_rna_count.csv'))
				
				if self.mode == 'translation':
					data_ordered_filtered = []
					order_pep_f = []
					peptides_order = []

					for peptide_type, peptides in self.peptides_by_type.items():
						for peptide in peptides:
							try:
								data_ordered_filtered.extend(data_filtered[peptide])
								peptides_order.append(peptide)

								if peptide_type not in order_pep_f:
									order_pep_f.append(peptide_type)
									order_f.append(1)
								else:
									ind = order_pep_f.index(peptide_type)
									order_f[ind] += 1
							except KeyError:
								pass

					df_counts_filtered = pd.DataFrame(data_ordered_filtered, columns=['Peptides', 'Alignments', 'BAM Files', 'Read Counts'])
					df_counts_filtered = df_counts_filtered.groupby(['Peptides', 'BAM Files'])['Read Counts'].sum().reset_index()
					df_counts_filtered = df_counts_filtered.pivot("Peptides", "BAM Files", "Read Counts").reset_index()
					sorterIndex = dict(zip(peptides_order, range(len(peptides_order))))
					df_counts_filtered['order'] = df_counts_filtered['Peptides'].map(sorterIndex)
					df_counts_filtered.sort_values(['order'], ascending = [True], inplace = True)
					df_counts_filtered.drop('order', 1, inplace = True)
					df_counts_filtered.set_index('Peptides',inplace=True)
					_thread.start_new_thread(self.save_info_counts, (df_counts_filtered, '', '_rna_ribo_count.csv'))
			
			t_2 = time.time()
			total = t_2-t_0
			logging.info('Total time run function get_counts to end : %f min', (total/60.0))
		else:
			logging.info('Count information already collected in the output folder : %s --> Skipping this step!', self.path_to_output_temps_folder+self.name_exp+'_rna_count.csv')
			df_counts = pd.read_csv(self.path_to_output_temps_folder+self.name_exp+'_rna_count.csv', index_col=0)
			
			if self.light:
				name_path = self.path_to_output_folder_alignments+'Alignments_information_light.dic'
			else:
				name_path = self.path_to_output_folder_alignments+'Alignments_information.dic'

			with open(name_path, 'rb') as fp:
				perfect_alignments = pickle.load(fp)

			if self.mode == 'translation':
				df_counts_filtered = pd.read_csv(self.path_to_output_temps_folder+self.name_exp+'_rna_ribo_count.csv', index_col=0)

			df_alignments = pd.read_csv(self.path_to_output_temps_folder+self.name_exp+'_rna_count_All_alignments.csv', index_col=0)

		return df_counts, perfect_alignments, df_counts_filtered, order, order_f, df_alignments


	def get_counts_sample(self, bam, peptide_alignment, sequences):
		
		name_sample = bam[0]
		bam_file = bam[1][0]
		library = bam[1][1]
		sens = bam[1][2]

		count = 0
		peptide = peptide_alignment.split('_')[0]
		alignment = peptide_alignment.split('_')[1]
		strand = peptide_alignment.split('_')[2]

		chr = alignment.split(':')[0]
		
		region_to_query = chr+':'+alignment.split(':')[1].split('-')[0]+'-'+alignment.split(':')[1].split('-')[-1]

		counts = self.get_depth_with_view(region_to_query, bam_file, library, sens, strand, sequences)

		to_return = [[peptide, alignment, name_sample, strand]]

		for index, count in enumerate(counts):
			sequence = sequences[index]
			to_return.append([count, sequence])
		 
		return to_return


	def save_info_counts(self, df, to_write, type_save):
		if to_write != '':
			with open(self.path_to_output_temps_folder+self.name_exp+type_save.split('.')[0]+'_All_alignments.csv', 'w') as csvFile:
				writer = csv.writer(csvFile)
				writer.writerows(to_write)

			csvFile.close()

		df.to_csv(self.path_to_output_temps_folder+self.name_exp+type_save, index=True, header=True)
		logging.info('Counts Information saved to : %s ', self.path_to_output_temps_folder+self.name_exp+type_save)


	def get_depth_with_view(self, region_to_query, bam_file, library, sens, strand, sequences):

		contReads_to_return = []

		library = library.lower()
		sens = sens.lower()
		
		if library == 'unstranded':
			try:
				count_1 = pysam.view("-F0X100", bam_file, region_to_query)
			except pysam.utils.SamtoolsError: 
				return -1

			for seq in sequences:
				contReads = 0
				rcmcs = uf.reverseComplement(seq)
				contReads += count_1.count(seq)
				contReads += count_1.count(rcmcs)
				contReads_to_return.append(contReads)

		elif library == 'single-end':

			if ((strand == '+' and sens == 'forward') or (strand == '-' and sens == 'reverse')):
				try:
					count_1 = pysam.view("-F0X110", bam_file, region_to_query)
				except pysam.utils.SamtoolsError: 
					return -1
			elif ((strand == '-' and sens == 'forward') or (strand == '+' and sens == 'reverse')):
				try:
					count_1 = pysam.view("-F0X100", "-f0X10", bam_file, region_to_query)
				except pysam.utils.SamtoolsError: 
					return -1

			for seq in sequences:
				contReads = 0
				contReads += count_1.count(seq)
				contReads_to_return.append(contReads)
				
		elif library == 'pair-end':
			if ((strand == '+' and sens == 'forward') or (strand == '-' and sens == 'reverse')):
				try:
					count_1 = pysam.view("-F0X100", "-f0X60", bam_file, region_to_query)
					count_2 = pysam.view("-F0X100", "-f0X90", bam_file, region_to_query)
				except pysam.utils.SamtoolsError: 
					return -1

			elif ((strand == '-' and sens == 'forward') or(strand == '+' and sens == 'reverse')):
				try:
					count_1 = pysam.view("-F0X100", "-f0X50", bam_file, region_to_query)
					count_2 = pysam.view("-F0X100", "-f0XA0", bam_file, region_to_query)
				except pysam.utils.SamtoolsError: 
					return -1

			count_1_split = count_1.split('\n')
			count_2_split = count_2.split('\n')
			
			for seq in sequences:
				contReads = 0
				reads_name = set()
				rcmcs = uf.reverseComplement(seq)

				for read in count_1_split:
					if seq in read :
						name = read.split('\t')[0]
						reads_name.add(name)

				for read in count_2_split:
					if rcmcs in read :
						name = read.split('\t')[0]
						reads_name.add(name)

				contReads = len(reads_name)
				contReads_to_return.append(contReads)

			# reading htslib https://github.com/DecodeGenetics/graphtyper/issues/57
		return contReads_to_return

	
