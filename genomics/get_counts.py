import warnings
warnings.filterwarnings("ignore")
import os, logging, time, subprocess, pickle, multiprocessing, os, _thread, csv, collections, pysam
import genomics.get_alignments as get_alig
import pandas as pd
from pathos.multiprocessing import ProcessPool
import utils.useful_functions as uf

__author__ = "Maria Virginia Ruiz Cuevas"
__email__ = "maria.virginia.ruiz.cuevas@umontreal.ca"


NUM_WORKERS =  multiprocessing.cpu_count()

path_to_lib = '/'.join(os.path.abspath(__file__).split('/')[:-3])+'/lib/'
genomePath = '/share/tsa_project/gtex_cram/cram_cache/Homo_sapiens.GRCh38.fa'

class GetCounts:

	def __init__(self, path_to_output_folder, name_exp, mode):
		self.path_to_output_temps_folder = path_to_output_folder+'res/temps_files/'
		self.name_exp = name_exp
		self.mode = mode
		self.path_to_output_folder_alignments = path_to_output_folder+'alignments/'


	def ribo_counts(self, perfect_alignments, bam_files_list):

		df_counts = pd.DataFrame({'A' : []})
		exists = os.path.exists(self.path_to_output_temps_folder+self.name_exp+'_ribo_count.csv')
		
		if not exists:
			t_0 = time.time()
			logging.info('Alignments on ribosome profiling information ')

			to_write = {} 
			
			perfect_alignments_to_return = {} 

			keys = []
			values = []
			data = []
			info_bams = []
			bams = []

			for name_sample, info_bam in bam_files_list.items():
				info_bams.append((name_sample,info_bam))
			
			def get_index(sample):
				index = list(bam_files_list.keys()).index(sample)
				return index + 1

			total_samples = len(bam_files_list.keys())

			keys = perfect_alignments.keys()
			values = perfect_alignments.values()
			pool = ProcessPool(nodes=NUM_WORKERS)

			bams.extend([info_bams]*len(keys))
			results = pool.map(self.get_counts_sample, bams, keys, values)

			for res in results:
				
				for count_align in res: 
					peptide = count_align[0]
					alignment = count_align[1]
					sample = count_align[2]
					count = count_align[3]
					strand = count_align[4]
					index = get_index(sample)

					try:
						peptide_info_to_write = to_write[peptide]

						try:
							peptide_info_to_write[alignment][index] = count

						except KeyError:
							counts = [0]*total_samples
							to_add = [strand]
							to_add.extend(counts)
							to_add[index] = count
							peptide_info_to_write[alignment] = to_add

					except KeyError:
						counts = [0]*total_samples
						to_add = [strand]
						to_add.extend(counts)
						to_add[index] = count
						to_write[peptide] = {alignment: to_add}

					data.append(count_align[:-2])

					key = peptide+'_'+alignment+'_'+sequence
					if len(perfect_alignments[key][-1]) == 0:
						perfect_alignments[key][-1] = [0]*total_samples

					perfect_alignments[key][-1][index-1] = count
					perfect_alignments_to_return[key] = perfect_alignments[key]


			header = ['Peptide', 'Position', 'Strand']
			header.extend(list(bam_files_list.keys()))
			
			to_write_list = [header]
			for peptide, info_peptide in to_write.items():
				to_add = []
				for alignment, info_counts in info_peptide.items():
					#if sum(info_counts[1:]) > 0:
					to_add =  [peptide, alignment]
					to_add.extend(info_counts)
					to_write_list.append(to_add)

			df_alignments = pd.DataFrame(to_write_list[1:], columns=to_write_list[0])

			if len(data) > 0 :
				df_counts = pd.DataFrame(data, columns=['Peptides', 'Alignments', 'BAM Files', 'Read Counts'])
				df_counts = df_counts.groupby(['Peptides', 'BAM Files'])['Read Counts'].sum().reset_index()
				df_counts = df_counts.pivot("Peptides", "BAM Files", "Read Counts")
				_thread.start_new_thread(self.save_info_counts, (df_counts, to_write_list, '_ribo_count.csv'))
			
			name_path = self.path_to_output_folder_alignments+'Alignments_ribo_information.dic'
			with open(name_path, 'wb') as handle:
				pickle.dump(perfect_alignments_to_return, handle, protocol=pickle.HIGHEST_PROTOCOL)

			name_path = self.path_to_output_folder_alignments+'Alignments_information.dic'
			with open(name_path, 'wb') as handle:
				pickle.dump(perfect_alignments, handle, protocol=pickle.HIGHEST_PROTOCOL)

			t_2 = time.time()
			total = t_2-t_0
			logging.info('Total time run function get_counts to end : %f min', (total/60.0))
		else:
			logging.info('Count information already collected in the output folder : %s --> Skipping this step!', self.path_to_output_temps_folder+self.name_exp+'_ribo_count.csv')
			
			df_counts = pd.read_csv(self.path_to_output_temps_folder+self.name_exp+'_ribo_count.csv', index_col=0)

			with open(self.path_to_output_folder_alignments+'Alignments_information.dic', 'rb') as fp:
				perfect_alignments = pickle.load(fp)
			
			df_alignments = pd.read_csv(self.path_to_output_temps_folder+self.name_exp+'_ribo_count'+'_All_alignments.csv', index_col=0)

		return perfect_alignments, df_counts, df_alignments


	def get_counts(self, perfect_alignments, bam_files_list):

		df_counts = pd.DataFrame()
		df_counts_filtered = pd.DataFrame()
		exist_rna = os.path.exists(self.path_to_output_temps_folder+self.name_exp+'_rna_count.csv')
		exist_ribo = True

		if self.mode == 'translation':
			with open(self.path_to_output_folder_alignments+'Alignments_ribo_information.dic', 'rb') as fp:
				perfect_alignments_ribo = pickle.load(fp)
			data_filtered = []
			exist_ribo = os.path.exists(self.path_to_output_temps_folder+self.name_exp+'_rna_ribo_count.csv')

		
		if not exist_rna and exist_ribo:
			t_0 = time.time()
			to_write = {} 
			
			keys = []
			values = []

			#for k in sorted(perfect_alignments, key=perfect_alignments.get, reverse=True):
			#	keys.append(k)
			#	values.append(perfect_alignments[k])

			data = []
			
			info_bams = []
			bams = []
			
			for name_sample, info_bam in bam_files_list.items():
				info_bams.append((name_sample,info_bam))
			
			def get_index(sample):
				index = list(bam_files_list.keys()).index(sample)
				return index + 1

			total_samples = len(bam_files_list.keys())

			keys = perfect_alignments.keys()
			values = perfect_alignments.values()
			
			pool = ProcessPool(nodes=NUM_WORKERS)

			bams.extend([info_bams]*len(keys))

			results = pool.map(self.get_counts_sample, bams, keys, values)

			for res in results:
				
				for count_align in res:
					peptide = count_align[0]
					alignment = count_align[1]
					sample = count_align[2]
					count = count_align[3]
					strand = count_align[4]
					sequence = count_align[5]

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

					data.append(count_align[:-2])

					key = peptide+'_'+alignment+'_'+sequence
					if len(perfect_alignments[key][-2]) == 0:
						perfect_alignments[key][-2] = [0]*total_samples

					perfect_alignments[key][-2][index-1] = count

					if self.mode == 'translation':
						try:
							count_ribo = perfect_alignments_ribo[key][-1][index-1]
							if count_ribo > 0:
								data_filtered.append(count_align[:-2])
						except KeyError:
							pass

			header = ['Peptide', 'Position', 'MCS', 'Strand']
			header.extend(bam_files_list.keys())
			to_write_list = [header]
			
			for peptide, info_peptide in to_write.items():
				to_add = []
				for alignment, info_counts in info_peptide.items():
					position = alignment.split('_')[0]
					MCS = alignment.split('_')[1]
					#if sum(info_counts[1:]) > 0:
					to_add =  [peptide, position, MCS]
					to_add.extend(info_counts)
					to_write_list.append(to_add)

			df_alignments = pd.DataFrame(to_write_list[1:], columns=to_write_list[0])

			# Replace dic
			name_path = self.path_to_output_folder_alignments+'/Alignments_information.dic'
			with open(name_path, 'wb') as handle:
				pickle.dump(perfect_alignments, handle, protocol=pickle.HIGHEST_PROTOCOL)


			if len(data) > 0 :
				df_counts = pd.DataFrame(data, columns=['Peptides', 'Alignments', 'BAM Files', 'Read Counts'])
				df_counts = df_counts.groupby(['Peptides', 'BAM Files'])['Read Counts'].sum().reset_index()
				df_counts = df_counts.pivot("Peptides", "BAM Files", "Read Counts")
				_thread.start_new_thread(self.save_info_counts, (df_counts, to_write_list, '_rna_count.csv'))
				
				if self.mode == 'translation':
					df_counts_filtered = pd.DataFrame(data_filtered, columns=['Peptides', 'Alignments', 'BAM Files', 'Read Counts'])
					df_counts_filtered = df_counts_filtered.groupby(['Peptides', 'BAM Files'])['Read Counts'].sum().reset_index()
					df_counts_filtered = df_counts_filtered.pivot("Peptides", "BAM Files", "Read Counts")
					_thread.start_new_thread(self.save_info_counts, (df_counts_filtered, '', '_rna_ribo_count.csv'))
			
			t_2 = time.time()
			total = t_2-t_0
			logging.info('Total time run function get_counts to end : %f min', (total/60.0))
		else:
			logging.info('Count information already collected in the output folder : %s --> Skipping this step!', self.path_to_output_temps_folder+self.name_exp+'_rna_count.csv')
			df_counts = pd.read_csv(self.path_to_output_temps_folder+self.name_exp+'_rna_count.csv', index_col=0)
			
			with open(self.path_to_output_folder_alignments+'/Alignments_information.dic', 'rb') as fp:
				perfect_alignments = pickle.load(fp)

			if self.mode == 'translation':
				df_counts_filtered = pd.read_csv(self.path_to_output_temps_folder+self.name_exp+'_rna_ribo_count.csv', index_col=0)

			df_alignments = pd.read_csv(self.path_to_output_temps_folder+self.name_exp+'_rna_count'+'_All_alignments.csv', index_col=0)

		return df_counts, perfect_alignments, df_counts_filtered, df_alignments


	def get_counts_sample(self, bams, peptide_alignment, info_alignment):
		to_return = []
		
		for info_bam in bams:
			count = 0
			name_sample = info_bam[0]
			bam_file = info_bam[1][0]
			library = info_bam[1][1]
			sens = info_bam[1][2]

			peptide = peptide_alignment.split('_')[0]
			alignment = peptide_alignment.split('_')[1]
			sequence = peptide_alignment.split('_')[2]

			for mcs in info_alignment[0]:
				strand = mcs[0]
				
				chr = alignment.split(':')[0]
				
				rcmcs = uf.reverseComplement(sequence)

				region_to_query = chr+':'+alignment.split(':')[1].split('-')[0]+'-'+alignment.split(':')[1].split('-')[-1]
				count = self.get_depth_with_view(region_to_query, bam_file, library, sens, strand, sequence)
				to_return.append([peptide, alignment, name_sample, count, strand, sequence])
		 
		return to_return

	def save_info_counts(self, df, to_write, type_save):
		if to_write != '':
			with open(self.path_to_output_temps_folder+self.name_exp+type_save.split('.')[0]+'_All_alignments.csv', 'w') as csvFile:
				writer = csv.writer(csvFile)
				writer.writerows(to_write)

			csvFile.close()

		df.to_csv(self.path_to_output_temps_folder+self.name_exp+type_save, index=True, header=True)
		logging.info('Counts Information saved to : %s ', self.path_to_output_temps_folder+self.name_exp+type_save)


	def get_depth_with_view(self, region_to_query, bam_file, library, sens, strand, seq):

		contReads = 0
		library = library.lower()
		sens = sens.lower()
		
		reads_name = set()
		rcmcs = uf.reverseComplement(seq)

		#if '.bam' in bam_file:
		if library == 'unstranded':

			count_1 = pysam.view("-F0X100", bam_file, region_to_query)
			contReads += count_1.count(seq)
			contReads += count_1.count(rcmcs)

		elif library == 'single-end':

			if ((strand == '+' and sens == 'forward') or (strand == '-' and sens == 'reverse')):
				count_1 = pysam.view("-F0X110", bam_file, region_to_query)
			elif ((strand == '-' and sens == 'forward') or (strand == '+' and sens == 'reverse')):
				count_1 = pysam.view("-F0X100", "-f0X10", bam_file, region_to_query)

			contReads += count_1.count(seq)

		elif library == 'pair-end':
			if ((strand == '+' and sens == 'forward') or (strand == '-' and sens == 'reverse')):
				count_1 = pysam.view("-F0X100", "-f0X60", bam_file, region_to_query)
				count_2 = pysam.view("-F0X100", "-f0X90", bam_file, region_to_query)

			elif ((strand == '-' and sens == 'forward') or(strand == '+' and sens == 'reverse')):
				count_1 = pysam.view("-F0X100", "-f0X50", bam_file, region_to_query)
				count_2 = pysam.view("-F0X100", "-f0XA0", bam_file, region_to_query)

			count_1_split = count_1.split('\n')
			count_2_split = count_2.split('\n')
			
			for read in count_1_split:
				if seq in read :
					name = read.split('\t')[0]
					reads_name.add(name)
			for read in count_2_split:
				if rcmcs in read :
					name = read.split('\t')[0]
					reads_name.add(name)

			contReads = len(reads_name)
			# reading htslib https://github.com/DecodeGenetics/graphtyper/issues/57
		return contReads

	
