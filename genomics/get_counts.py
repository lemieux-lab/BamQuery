import warnings, os
warnings.filterwarnings("ignore")
import os, time, subprocess, pickle, multiprocessing, os, _thread, csv, collections, pysam, copy
import genomics.get_counts_from_sample as get_counts_sample
import pandas as pd
from pathos.multiprocessing import ProcessPool
import utils.useful_functions as uf
import numpy as np

__author__ = "Maria Virginia Ruiz Cuevas"
__email__ = "maria.virginia.ruiz.cuevas@umontreal.ca"


NUM_WORKERS =  multiprocessing.cpu_count()

path_to_lib = '/'.join(os.path.abspath(__file__).split('/')[:-3])+'/lib/'

class GetCounts:

	def __init__(self, path_to_output_folder, name_exp, mode, light, peptides_by_type, super_logger):
		self.light = light
		if self.light:
			self.path_to_output_temps_folder = path_to_output_folder+'res_light/temps_files/'
		else:
			self.path_to_output_temps_folder = path_to_output_folder+'res/temps_files/'
		
		self.name_exp = name_exp
		self.mode = mode
		self.path_to_output_folder_alignments = path_to_output_folder+'alignments/'
		self.peptides_by_type = peptides_by_type
		self.super_logger = super_logger
		

	def ribo_counts(self, perfect_alignments, bam_files_list):

		df_counts = pd.DataFrame()
		exists = os.path.exists(self.path_to_output_temps_folder+self.name_exp+'_ribo_count.csv')
		
		order = []

		last_treated_bam_file = os.path.exists(self.path_to_output_folder_alignments+'info_trated_ribo_bam_files.pkl')

		if last_treated_bam_file:

			with open(self.path_to_output_folder_alignments+'info_trated_ribo_bam_files.pkl', 'rb') as fp:
				last_treated_bam_file = pickle.load(fp)
		else:
			last_treated_bam_file = 0

		if self.light:
			alignment_information = self.path_to_output_folder_alignments+'Alignments_information_light.dic'
			alignment_information_ribo = self.path_to_output_folder_alignments+'Alignments_ribo_information_light.dic'
		else:
			alignment_information = self.path_to_output_folder_alignments+'Alignments_information.dic'
			alignment_information_ribo = self.path_to_output_folder_alignments+'Alignments_ribo_information.dic'

		times = []

		if not exists:
			t_0 = time.time()
			self.super_logger.info('Alignments on ribosome profiling information ')

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
			
			self.super_logger.info('Total MCS mapped : %s ', str(len(keys)))

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

			self.super_logger.info('Total unique regions : %s ', str(len(keys)))
			
			pool = ProcessPool(nodes=NUM_WORKERS)

			for idx, bam_file in enumerate(info_bams):

				if idx > last_treated_bam_file:

					t0_bam_file = time.time()
					bams = [bam_file] * len(keys)
					results = pool.map(get_counts_sample.get_counts_sample, bams, keys, values)
					not_permission = False

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
								if count == -1:
									not_permission = True
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
					times.append(time_final)

					if not_permission:
						self.super_logger.info('Bam File : %s %s couldn\'t be processed. Failed to open, permission denied. Time : %f min', str(idx), bam_file[0], time_final)
					else:
						self.super_logger.info('Processed Bam File : %s %s. Time : %f min', str(idx), bam_file[0], time_final)

					if (idx % 50) == 0:
						with open(self.path_to_output_folder_alignments+'info_trated_ribo_bam_files.pkl', 'wb') as f:  
							pickle.dump(idx, f)

						with open(alignment_information, 'wb') as handle:
							pickle.dump(perfect_alignments, handle, protocol=pickle.HIGHEST_PROTOCOL)

						with open(alignment_information_ribo, 'wb') as handle:
							pickle.dump(perfect_alignments_to_return, handle, protocol=pickle.HIGHEST_PROTOCOL)
				else:
					self.super_logger.info('Bam File already processed: %s %s.', str(idx), bam_file[0])

			pool.close()
			pool.join()
			pool.clear()
			
			with open(self.path_to_output_folder_alignments+'info_trated_ribo_bam_files.pkl', 'wb') as f:  
				pickle.dump(idx, f)

			with open(alignment_information, 'wb') as handle:
				pickle.dump(perfect_alignments, handle, protocol=pickle.HIGHEST_PROTOCOL)

			with open(alignment_information_ribo, 'wb') as handle:
				pickle.dump(perfect_alignments_to_return, handle, protocol=pickle.HIGHEST_PROTOCOL)

			self.super_logger.info('Average time to process a BamFile : %f min', np.mean(times))

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
			self.super_logger.info('Total time run function get_counts to end : %f min', (total/60.0))
		else:
			self.super_logger.info('Count information already collected in the output folder : %s --> Skipping this step!', self.path_to_output_temps_folder+self.name_exp+'_ribo_count.csv')
			df_counts = pd.read_csv(self.path_to_output_temps_folder+self.name_exp+'_ribo_count.csv', index_col=0)

			with open(alignment_information, 'rb') as fp:
				perfect_alignments = pickle.load(fp)
			
			df_alignments = pd.read_csv(self.path_to_output_temps_folder+self.name_exp+'_ribo_count_All_alignments.csv', index_col=0)

		return perfect_alignments, df_counts, order, df_alignments


	def get_counts(self, perfect_alignments, bam_files_list):

		df_counts = pd.DataFrame()
		df_counts_filtered = pd.DataFrame()

		exist_rna = os.path.exists(self.path_to_output_temps_folder+self.name_exp+'_rna_count.csv')
		exist_ribo = True

		last_treated_bam_file = os.path.exists(self.path_to_output_folder_alignments+'info_trated_bam_files.pkl')

		if last_treated_bam_file:

			with open(self.path_to_output_folder_alignments+'info_trated_bam_files.pkl', 'rb') as fp:
				last_treated_bam_file = pickle.load(fp)

			try:
				with open(self.path_to_output_folder_alignments+'data.pkl', 'rb') as fp:
					data = pickle.load(fp)

				with open(self.path_to_output_folder_alignments+'to_write.pkl', 'rb') as fp:
					to_write = pickle.load(fp)
			except FileNotFoundError:
				to_write = {} 
				data = {}

		else:
			last_treated_bam_file = -1
			to_write = {} 
			data = {}

		if self.light:
			alignment_information = self.path_to_output_folder_alignments+'Alignments_information_light.dic'
		else:
			alignment_information = self.path_to_output_folder_alignments+'Alignments_information.dic'
			
		if self.mode == 'translation':

			if self.light:
				name_path = self.path_to_output_folder_alignments+'Alignments_ribo_information.dic'
			else:
				name_path = self.path_to_output_folder_alignments+'Alignments_ribo_information_light.dic'

			with open(name_path, 'rb') as fp:
				perfect_alignments_ribo = pickle.load(fp)

			exist_ribo = os.path.exists(self.path_to_output_temps_folder+self.name_exp+'_rna_ribo_count.csv')
			try:
				with open(self.path_to_output_folder_alignments+'data_filtered.pkl', 'rb') as fp:
					data_filtered = pickle.load(fp)
			except:
				data_filtered = {}
			

		order_f = []
		order = []
		times = []

		if (not exist_rna and exist_ribo):
			t_0 = time.time()
			
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
			self.super_logger.info('Total MCS mapped : %s ', str(len(keys)))

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
			
			self.super_logger.info('Total unique regions : %s ', str(len(keys)))

			pool = ProcessPool(nodes=NUM_WORKERS)

			for idx, bam_file in enumerate(info_bams):
				
				if idx > last_treated_bam_file:

					t0_bam_file = time.time()
					bams = [bam_file] * len(keys)
					results = pool.map(get_counts_sample.get_counts_sample, bams, keys, values)
					not_permission = False

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
								if count == -1:
									not_permission = True
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
					times.append(time_final)

					if not_permission:
						self.super_logger.info('Bam File : %s %s couldn\'t be processed. Failed to open, permission denied. Time : %f min', str(idx), bam_file[0], time_final)
					else:
						self.super_logger.info('Processed Bam File : %s %s. Time : %f min', str(idx), bam_file[0], time_final)

					if (idx % 100 == 0) and (idx != 0):
						self.super_logger.info('Saving information for Bam Files processed')

						with open(self.path_to_output_folder_alignments+'info_trated_bam_files.pkl', 'wb') as f:  
							pickle.dump(idx, f)

						with open(alignment_information, 'wb') as handle:
							pickle.dump(perfect_alignments, handle, protocol=pickle.HIGHEST_PROTOCOL)

						with open(self.path_to_output_folder_alignments+'to_write.pkl', 'wb') as handle:
							pickle.dump(to_write, handle, protocol=pickle.HIGHEST_PROTOCOL)

						with open(self.path_to_output_folder_alignments+'data.pkl', 'wb') as handle:
							pickle.dump(data, handle, protocol=pickle.HIGHEST_PROTOCOL)

						if self.mode == 'translation':
							with open(self.path_to_output_folder_alignments+'data_filtered.pkl', 'wb') as handle:
								pickle.dump(data_filtered, handle, protocol=pickle.HIGHEST_PROTOCOL)							
				else:
					self.super_logger.info('Bam File already processed: %s %s.', str(idx), bam_file[0])

			pool.close()
			pool.join()
			pool.clear()

			self.super_logger.info('Average time to process a BamFile : %f min', np.mean(times))

			with open(self.path_to_output_folder_alignments+'info_trated_bam_files.pkl', 'wb') as f:  
				pickle.dump(idx, f)
			
			with open(alignment_information, 'wb') as handle:
				pickle.dump(perfect_alignments, handle, protocol=pickle.HIGHEST_PROTOCOL)


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

			order_pep = []
			peptides_order = []
			
			for peptide_type, peptides in self.peptides_by_type.items():
				for peptide in peptides:
					peptides_order.append(peptide)
					if peptide_type not in order_pep:
						order_pep.append(peptide_type)
						order.append(1)
					else:
						ind = order_pep.index(peptide_type)
						order[ind] += 1
			
			# https://stackoverflow.com/questions/23482668/sorting-by-a-custom-list-in-pandas

			df_counts = df_alignments.groupby(['Peptide']).sum().reset_index()

			sorterIndex = dict(zip(peptides_order, range(len(peptides_order))))
			df_counts['order'] = df_counts['Peptide'].map(sorterIndex)
			df_counts.sort_values(['order'], ascending = [True], inplace = True)
			df_counts.drop('order', 1, inplace = True)
			df_counts.set_index('Peptide',inplace=True)

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

				os.remove(self.path_to_output_folder_alignments+'data_filtered.pkl')
			
			t_2 = time.time()
			total = t_2-t_0
			self.super_logger.info('Total time run function get_counts to end : %f min', (total/60.0))
		else:
			self.super_logger.info('Count information already collected in the output folder : %s --> Skipping this step!', self.path_to_output_temps_folder+self.name_exp+'_rna_count.csv')
			df_counts = pd.read_csv(self.path_to_output_temps_folder+self.name_exp+'_rna_count.csv', index_col=0)
			
			with open(alignment_information, 'rb') as fp:
				perfect_alignments = pickle.load(fp)

			if self.mode == 'translation':
				df_counts_filtered = pd.read_csv(self.path_to_output_temps_folder+self.name_exp+'_rna_ribo_count.csv', index_col=0)

			df_alignments = pd.read_csv(self.path_to_output_temps_folder+self.name_exp+'_rna_count_All_alignments.csv', index_col=0)

		try:
			os.remove(self.path_to_output_folder_alignments+'to_write.pkl')
			os.remove(self.path_to_output_folder_alignments+'data.pkl')
		except FileNotFoundError:
			pass
		
		return df_counts, perfect_alignments, df_counts_filtered, order, order_f, df_alignments


	def save_info_counts(self, df, to_write, type_save):
		if to_write != '':
			with open(self.path_to_output_temps_folder+self.name_exp+type_save.split('.')[0]+'_All_alignments.csv', 'w') as csvFile:
				writer = csv.writer(csvFile)
				writer.writerows(to_write)

			csvFile.close()

		df.to_csv(self.path_to_output_temps_folder+self.name_exp+type_save, index=True, header=True)
		self.super_logger.info('Counts Information saved to : %s ', self.path_to_output_temps_folder+self.name_exp+type_save)



	
