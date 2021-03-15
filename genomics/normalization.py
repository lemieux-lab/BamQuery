import os, logging, time, pickle, multiprocessing, _thread, csv, math
import pandas as pd

NUM_WORKERS =  multiprocessing.cpu_count()

__author__ = "Maria Virginia Ruiz Cuevas"

class Normalization:

	def __init__(self, path_to_output_folder, name_exp, bam_files_list, mode):
		path_to_lib = '/'.join(os.path.abspath(__file__).split('/')[:-3])+'/lib/'
		self.path_to_output_folder = path_to_output_folder+'res/'
		self.name_exp = name_exp
		self.path_to_all_counts_file = path_to_lib+"allcounts.dic"
		self.bam_files_list = bam_files_list
		self.mode = mode

	def get_normalization(self, df_counts, type_save):

		exists = os.path.exists(self.path_to_output_folder+self.name_exp+type_save)
		
		if not exists:
		
			t_0 = time.time()
			count_reads = {}

			indexes_column_names = ['Peptides']
			indexes_column_names.extend(list(df_counts.columns))
			norm_matrix = []

			try:
				with open(self.path_to_all_counts_file, 'rb') as fp:
					dictionary_total_reads_bam_files = pickle.load(fp)
			except FileNotFoundError:
				print ('Allcounts dictionary is not found in the lib path. If this is the first time you are running \
					BAMQuery or you are querying new samples, this will take a bit of time. Go for coffee! If it is not \
					the case you should be worried, the Allcounts dictionary have been lost, now is generating a new Allcounts dictionary \
					with the samples in this query.')
				logging.info('Allcounts dictionary is not found in the lib path. If this is the first time you are running \
					BAMQuery or you are querying new samples, this will take a bit of time. Go for coffee! If it is not \
					the case you should be worried, the Allcounts dictionary have been lost, now is generating a new Allcounts dictionary \
					with the samples in this query.')

			while len(count_reads) != len(self.bam_files_list):
				try:
					with open(self.path_to_all_counts_file, 'rb') as fp:
						dictionary_total_reads_bam_files = pickle.load(fp)

					for bam_file, bamfile_info in self.bam_files_list.items():
						try:
							info_bam_file = dictionary_total_reads_bam_files[bam_file]
							count_reads_saved = info_bam_file[1]
							count_reads[bam_file] = count_reads_saved
						except KeyError:
							pass
				except FileNotFoundError:
					pass
				
			for index, row in df_counts.iterrows():
				aux = [0.0]*len(indexes_column_names)
				peptide = row.name
				aux[0] =  peptide
				for bam_file, count_bam_file in count_reads.items():
					df_counts[bam_file].astype('float')
					count = int(row[bam_file])
					norm = math.log(((count/(int(count_bam_file)*1.0))*100000000.0)+1,10)
					row[bam_file] = norm
					index = indexes_column_names.index(bam_file)
					aux[index] = norm
				norm_matrix.append(aux)

			def_norm = pd.DataFrame(norm_matrix, columns = indexes_column_names)
			def_norm = def_norm.set_index(['Peptides'])

			_thread.start_new_thread(self.save_info_counts, (def_norm, type_save))
			
			t_2 = time.time()
			total = t_2-t_0
			logging.info('Total time run function get_normalization to end : %f min', (total/60.0))
		else:
			logging.info('Normalization information already collected in the output folder : %s --> Skipping this step!', self.path_to_output_folder+self.name_exp+type_save)
			def_norm = pd.read_csv(self.path_to_output_folder+self.name_exp+type_save, index_col=0)
		
		return def_norm

	def save_info_counts(self, df, type_save):
		df.to_csv(self.path_to_output_folder+self.name_exp+type_save, index=True, header=True)
		logging.info('Normalization Information saved to : %s ', self.path_to_output_folder+self.name_exp+type_save)

