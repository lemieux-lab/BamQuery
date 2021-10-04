import os, logging, time, pickle, subprocess, multiprocessing, _thread, csv, math
import pandas as pd
import numpy as np


NUM_WORKERS =  multiprocessing.cpu_count()

__author__ = "Maria Virginia Ruiz Cuevas"

class Normalization:

	def __init__(self, path_to_output_folder, name_exp, bam_files_list, peptides_types, mode, light):
		path_to_lib = '/'.join(os.path.abspath(__file__).split('/')[:-3])+'/lib/'
		self.light = light
		
		if self.light:
			self.path_to_output_folder = path_to_output_folder+'res_light/'
			self.path_to_output_temps_folder = path_to_output_folder+'res_light/temps_files/'
			self.path_to_output_aux_folder = path_to_output_folder+'res_light/AUX_files/'
			self.path_to_output_aux_processed_folder = path_to_output_folder+'res_light/AUX_files/processed/'
		else:
			self.path_to_output_folder = path_to_output_folder+'res/'
			self.path_to_output_temps_folder = path_to_output_folder+'res/temps_files/'
			self.path_to_output_aux_folder = path_to_output_folder+'res/AUX_files/'
			self.path_to_output_aux_processed_folder = path_to_output_folder+'res/AUX_files/processed/'
			

		self.name_exp = name_exp
		self.path_to_all_counts_file = path_to_lib+"allcounts.dic"
		self.path_to_tissues_file = path_to_lib+"tissues.dic"
		self.bam_files_list = bam_files_list
		self.peptides_types = peptides_types
		self.mode = mode


	def get_normalization(self, df_counts, type_save):

		exists = os.path.exists(self.path_to_output_temps_folder+self.name_exp+type_save)
		
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


			try:
				with open(self.path_to_tissues_file, 'rb') as fp:
					dictionary_tissues_bam_files = pickle.load(fp)
			except ValueError:
				import pickle5
				with open(self.path_to_tissues_file, 'rb') as fp:
					dictionary_tissues_bam_files = pickle5.load(fp)

			if os.path.exists(self.path_to_output_aux_folder+"bam_files_tissues.csv"):
				
				df = pd.read_csv(self.path_to_output_aux_folder+"bam_files_tissues.csv", header = 0)
				for index, row in df.iterrows():
					sample = row['Sample']
					tissue_name = row['Tissue']
					tissue_type = row['Tissue_type']
					short_list = row['Short_list'].lower()
					dictionary_tissues_bam_files[sample] = [tissue_name, tissue_type, short_list]
				
				
				with open(self.path_to_tissues_file, 'wb') as handle:
					pickle.dump(dictionary_tissues_bam_files, handle, protocol=pickle.HIGHEST_PROTOCOL)
				
			not_in = set()
			tissues_for_samples = {}

			while len(count_reads) != len(self.bam_files_list):
				try:
					try:
						with open(self.path_to_tissues_file, 'rb') as fp:
							dictionary_tissues_bam_files = pickle.load(fp)
					except ValueError:
						import pickle5
						with open(self.path_to_tissues_file, 'rb') as fp:
							dictionary_tissues_bam_files = pickle5.load(fp)
					
					for bam_file, bamfile_info in self.bam_files_list.items():

						try:
							info_bam_file = dictionary_total_reads_bam_files[bam_file]
							count_reads_saved = info_bam_file[1]
							count_reads[bam_file] = count_reads_saved
						except KeyError:
							pass

						try:
							tissue_info = dictionary_tissues_bam_files[bam_file]
							tissue = tissue_info[0]
							tissue_type = tissue_info[1]
							short_list = tissue_info[2]
							try:
								tissues_for_samples[tissue][0].append(bam_file)
							except KeyError:
								tissues_for_samples[tissue] = [[bam_file], tissue_type, short_list]

						except KeyError:
							not_in.add(bam_file)

				except :
					pass

			data = [['Sample category', 'sample_ids', 'Project', 'short_list']]
			for tissue, info_tissue in tissues_for_samples.items():
				aux = []
				aux.append(tissue)
				bam_files = ' '.join(info_tissue[0])
				aux.append(bam_files)
				aux.append(info_tissue[1])
				aux.append(info_tissue[2])
				data.append(aux)

			with open(self.path_to_output_folder+"info_bam_files_tissues.csv", 'w') as csvFile:
				writer = csv.writer(csvFile)
				writer.writerows(data)

			if len(not_in) > 0 :
				data = [['Sample', 'Tissue', 'Tissue_type', 'Short_list']]
				for sample in not_in:
					to_add = ['']*4
					to_add[0] = sample
					data.append(to_add)

				with open(self.path_to_output_aux_folder+"bam_files_tissues.csv", 'w') as csvFile:
					writer = csv.writer(csvFile)
					writer.writerows(data)

				raise Exception("\nBefore to continue you must provide the tissue type for the bam files annotated in the file : "+ self.path_to_output_aux_folder+"bam_files_tissues.csv. Please enter for each sample : tissue, tissue_type, shortlist." )
			
			info_tissue_peptide = {}
			
			for index, row in df_counts.iterrows():
				aux = [0.0]*len(indexes_column_names)
				peptide = row.name
				aux[0] =  peptide
				for bam_file, count_bam_file in count_reads.items():
					#df_counts[bam_file].astype('float')
					count = int(row[bam_file])
					norm = math.log(((count/(int(count_bam_file)*1.0))*100000000.0)+1,10)
					#row[bam_file] = norm
					index = indexes_column_names.index(bam_file)
					aux[index] = norm
					
					tissue_info = dictionary_tissues_bam_files[bam_file]
					tissue = tissue_info[0]
					
					try:
						info_tissue = info_tissue_peptide[tissue]
						try:
							info_tissue[peptide][0].append(norm)
						except KeyError:
							info_tissue[peptide] = [[norm]]
					except KeyError:
						info_tissue_peptide[tissue] = {peptide: [[norm]]}
					

				norm_matrix.append(aux)

			exp = type_save[1:].split('.csv')[0]+'/'
			for tissue, peptides in info_tissue_peptide.items():
				to_write = 'Peptide\tPeptide_type\tTissue\tTissue_type\tShort_list\tmedian\tmean\n'
				for peptide, info in peptides.items():
					mean = np.mean(info[0])
					median = np.median(info[0])
					info.append(mean)
					info.append(median)
					peptide_type = self.peptides_types[peptide]
					tissue_type = tissues_for_samples[tissue][1]
					short_list = tissues_for_samples[tissue][2]
					to_write += peptide+'\t'+peptide_type+'\t'+tissue+'\t'+tissue_type+'\t'+short_list+'\t'+str(median)+'\t'+str(mean)+'\n'

				with open(self.path_to_output_aux_processed_folder+exp+tissue+"_processed.txt", 'w') as f:
					f.write(to_write)
					f.close()

			def_norm = pd.DataFrame(norm_matrix, columns = indexes_column_names)
			def_norm = def_norm.set_index(['Peptides'])

			_thread.start_new_thread(self.save_info_counts, (def_norm, type_save))
			
			t_2 = time.time()
			total = t_2-t_0
			logging.info('Total time run function get_normalization to end : %f min', (total/60.0))
		else:
			logging.info('Normalization information already collected in the output folder : %s --> Skipping this step!', self.path_to_output_temps_folder+self.name_exp+type_save)
			def_norm = pd.read_csv(self.path_to_output_temps_folder+self.name_exp+type_save, index_col=0)
		
		return def_norm

	def save_info_counts(self, df, type_save):
		df.to_csv(self.path_to_output_temps_folder+self.name_exp+type_save, index=True, header=True)
		logging.info('Normalization Information saved to : %s ', self.path_to_output_temps_folder+self.name_exp+type_save)

# 				t1 				ts2 			ts3 			ts4
# rphm			x				x				x				x
# log(rphm)		log(x+1,10)		log(x+1,10)		log(x+1,10)		log(x+1,10)
# mean_log		sum(log(x+1,10)	+ log(x+1,10)+log(x+1,10)+log(x+1,10))/4
# arithmetic_mean	sum(x1	+ x2 + x3 + x4)/4

