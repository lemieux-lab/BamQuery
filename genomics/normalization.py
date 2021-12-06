import os, time, pickle, subprocess, multiprocessing, _thread, csv, math, copy
import pandas as pd
import numpy as np


NUM_WORKERS =  multiprocessing.cpu_count()

__author__ = "Maria Virginia Ruiz Cuevas"

class Normalization:

	def __init__(self, path_to_output_folder, name_exp, peptides_types, mode, light, super_logger):
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
		self.path_to_all_counts_file = path_to_lib+"Bam_files_info.dic"
		self.peptides_types = peptides_types
		self.mode = mode
		self.super_logger = super_logger


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
			except ValueError:
				import pickle5
				with open(self.path_to_all_counts_file, 'rb') as fp:
					dictionary_total_reads_bam_files = pickle5.load(fp)
			
			# 0: path, 1: number, 2: Tissue, 3: Tissue_type, 4: Shortlist, 5: sequencing, 6: library, 7: user
			if os.path.exists(self.path_to_output_aux_folder+"bam_files_tissues.csv"):
				self.super_logger.info('Adding tissue information to Bam Files !')
				
				df = pd.read_csv(self.path_to_output_aux_folder+"bam_files_tissues.csv", header = 0)
				for index, row in df.iterrows():
					try:
						sample = str(row['Sample'])
						tissue_name = str(row['Tissue'])
						tissue_type = str(row['Tissue_type'])
						short_list = str(row['Short_list']).lower()
						if sample != 'nan' and tissue_name != 'nan' and tissue_type != 'nan':
							dictionary_total_reads_bam_files[sample][2] = tissue_name
							dictionary_total_reads_bam_files[sample][3] = tissue_type
							dictionary_total_reads_bam_files[sample][4] = short_list
						else:
							raise Exception("\nBefore to continue you must provide the tissue type for the bam files annotated in the file : "+ self.path_to_output_aux_folder+"bam_files_tissues.csv. Please enter for each sample : tissue, tissue_type, shortlist." )
					except:
						raise Exception("\nBefore to continue you must provide the tissue type for the bam files annotated in the file : "+ self.path_to_output_aux_folder+"bam_files_tissues.csv. Please enter for each sample : tissue, tissue_type, shortlist." )
				with open(self.path_to_all_counts_file, 'wb') as handle:
					pickle.dump(dictionary_total_reads_bam_files, handle, protocol=pickle.HIGHEST_PROTOCOL)
				
			tissues_for_samples = {}
			bam_files_list = list(df_counts.columns)
			counts_bam_files = []
			
			for name_bam_file in bam_files_list:

				try:
					info_bam_file = dictionary_total_reads_bam_files[name_bam_file]
					path = info_bam_file[0]
					count = info_bam_file[1]
					tissue = info_bam_file[2]
					tissue_type = info_bam_file[3]
					shortlist = info_bam_file[4]
					sequencing = info_bam_file[5]
					library = info_bam_file[6]
					user = info_bam_file[7]
					counts_bam_files.append(count*1.0)

					if count == 0:
						raise Exception("\nBefore to continue you need to verify that the primary read count for the bam file "+name_bam_file+" is already included in the dictionary. To do so: verify in the log Get_Read_Count_BAM_directories.log that the primary read count processes have finished. Please re-launch BamQuery, once all the primary read counts have been included." )
					
					if tissue == '' or tissue_type == '' or shortlist == '':
						print (name_bam_file, info_bam_file)
						raise Exception("\nBefore to continue you must provide the tissue type for the bam files annotated in the file : "+ self.path_to_output_aux_folder+"bam_files_tissues.csv. Please enter for each sample : tissue, tissue_type, shortlist." )
					try:
						tissues_for_samples[tissue][0].append(name_bam_file)
					except KeyError:
						tissues_for_samples[tissue] = [[name_bam_file], tissue_type, shortlist]
				except KeyError:
					raise Exception("\nBefore to continue you need to verify that the primary read count for the bam file "+name_bam_file+" is already included in the dictionary. To do so: verify in the log Get_Read_Count_BAM_directories.log that the primary read count processes have finished. Please re-launch BamQuery, once all the primary read counts have been included." )
					
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

			info_tissue_peptide = {}
			
			def_norm = copy.deepcopy(df_counts)

			def normalize(row):
				new_row = ((row / counts_bam_files)*100000000.0)+1
				return new_row

			def_norm = def_norm.apply(lambda row : normalize(row), axis = 1)
			def_norm_log = def_norm.apply(np.log10, axis = 1)
			
			for index, row in def_norm_log.iterrows():
				peptide = row.name
				for bam_file_name, norm in row.items():
					tissue_name = dictionary_total_reads_bam_files[bam_file_name][2]
					tissue_type = dictionary_total_reads_bam_files[bam_file_name][3]
					short_list = dictionary_total_reads_bam_files[bam_file_name][4]

					if tissue_name == '':
						print (bam_file_name)
					key = (tissue_name, tissue_type, short_list)
					try:
						info_tissue = info_tissue_peptide[key]
						try:
							info_tissue[peptide][0].append(norm)
						except KeyError:
							info_tissue[peptide] = [[norm]]
					except KeyError:
						info_tissue_peptide[key] = {peptide: [[norm]]}

			exp = type_save[1:].split('.csv')[0]+'/'

			for key_info_tissue, peptides in info_tissue_peptide.items():
				
				to_write = 'Peptide\tPeptide_type\tTissue\tTissue_type\tShort_list\tmedian\tmean\n'
				tissue_name = key_info_tissue[0]
				tissue_type = key_info_tissue[1]
				short_list = key_info_tissue[2]

				for peptide, info in peptides.items():
					mean = np.mean(info[0])
					median = np.median(info[0])
					peptide_type = self.peptides_types[peptide]
					to_write += peptide+'\t'+peptide_type+'\t'+tissue_name+'\t'+tissue_type+'\t'+short_list+'\t'+str(median)+'\t'+str(mean)+'\n'

				with open(self.path_to_output_aux_processed_folder+exp+tissue_name+"_processed.txt", 'w') as f:
					f.write(to_write)
					f.close()
			_thread.start_new_thread(self.save_info_counts, (def_norm_log, type_save))
			
			t_2 = time.time()
			total = t_2-t_0
			self.super_logger.info('Total time run function get_normalization to end : %f min', (total/60.0))
		else:
			self.super_logger.info('Normalization information already collected in the output folder : %s --> Skipping this step!', self.path_to_output_temps_folder+self.name_exp+type_save)
			def_norm_log = pd.read_csv(self.path_to_output_temps_folder+self.name_exp+type_save, index_col=0)
		
		return def_norm_log

	def save_info_counts(self, df, type_save):
		df.to_csv(self.path_to_output_temps_folder+self.name_exp+type_save, index=True, header=True)
		self.super_logger.info('Normalization Information saved to : %s ', self.path_to_output_temps_folder+self.name_exp+type_save)

# 				t1 				ts2 			ts3 			ts4
# rphm			x				x				x				x
# log(rphm)		log(x+1,10)		log(x+1,10)		log(x+1,10)		log(x+1,10)
# mean_log		sum(log(x+1,10)	+ log(x+1,10)+log(x+1,10)+log(x+1,10))/4
# arithmetic_mean	sum(x1	+ x2 + x3 + x4)/4

