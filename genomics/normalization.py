import os, time, pickle, subprocess, multiprocessing, _thread, csv, math, copy
import pandas as pd
import numpy as np


NUM_WORKERS =  multiprocessing.cpu_count()

path_to_lib = '/'.join(os.path.abspath(__file__).split('/')[:-3])+'/lib/'
__author__ = "Maria Virginia Ruiz Cuevas"

class Normalization:

	def __init__(self, path_to_output_folder, name_exp, peptides_types, mode, light, super_logger, dev):
		
		self.light = light
		self.mode = mode
		
		if self.light:
			self.path_to_output_folder = path_to_output_folder+'res_light/'
			self.path_to_output_temps_folder = path_to_output_folder+'res_light/temps_files/'
		else:
			self.path_to_output_folder = path_to_output_folder+'res/'
			self.path_to_output_temps_folder = path_to_output_folder+'res/temps_files/'
			self.path_to_output_total_transcription_expression_heatmap_folder = path_to_output_folder+'plots/heat_maps/transcription_evidence_heatmap/average_transcription_expression_heatmap/'
		if self.mode == 'translation':
			self.path_to_output_folder = path_to_output_folder+'res_translation/'
			self.path_to_output_temps_folder = path_to_output_folder+'res_translation/temps_files/'
			self.path_to_output_total_transcription_expression_heatmap_folder = path_to_output_folder+'plots/heat_maps/translation_evidence_heatmap/average_translation_expression_heatmap/'
		
		self.name_exp = name_exp
		self.path_to_all_counts_file = path_to_lib+"Bam_files_info.dic"
		self.peptides_types = peptides_types
		self.super_logger = super_logger
		self.dev = dev

	def get_normalization(self, df_counts, type_save):

		exists = os.path.exists(self.path_to_output_temps_folder+self.name_exp+type_save)
		
		if not exists:
		
			t_0 = time.time()
			count_reads = {}

			norm_matrix = []

			incomplete = True

			while incomplete:

				try:
					with open(self.path_to_all_counts_file, 'rb') as fp:
						dictionary_total_reads_bam_files = pickle.load(fp)
				except ValueError:
					import pickle5
					with open(self.path_to_all_counts_file, 'rb') as fp:
						dictionary_total_reads_bam_files = pickle5.load(fp)
				
				tissues_for_samples = {}
				bam_files_list = list(df_counts.columns[2:])
				counts_bam_files = []
				
				all_counts_included = True

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

						if tissue == '' or tissue_type == '' or shortlist == '':
							raise Exception("\nBefore to continue you must provide the tissue type for the bam files annotated in the file : "+ self.path_to_output_aux_folder+"bam_files_tissues.csv. Please enter for each sample : tissue, tissue_type, shortlist." )
						try:
							tissues_for_samples[tissue][0].append(name_bam_file)
						except KeyError:
							tissues_for_samples[tissue] = [[name_bam_file], tissue_type, shortlist]

						if count == 0:
							all_counts_included = False
							break
						
					except KeyError:
						raise Exception("\nBefore to continue you need to verify that the primary read count for the bam file "+name_bam_file+" is already included in the dictionary. To do so: verify in the log Get_Read_Count_BAM_directories.log that the primary read count processes have finished. Please re-launch BamQuery, once all the primary read counts have been included." )
				
				if not all_counts_included:
					incomplete = True
					time.sleep(10)
				else:
					incomplete = False

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
			df_counts.set_index(['Peptide Type','Peptide'], inplace=True)
			def_norm = copy.deepcopy(df_counts)
			
			def normalize(row):
				row = ((row / counts_bam_files)*100000000.0)+1
				return row

			def_norm = def_norm.apply(lambda row : normalize(row), axis = 1)
			def_norm_log = def_norm.apply(np.log10, axis = 1)
			
			if not self.light:
				data = []
				for index, row in def_norm_log.iterrows():
					peptide_type = row.name[0]
					peptide = row.name[1]
					
					for bam_file_name, norm in row.items():
						tissue_name = dictionary_total_reads_bam_files[bam_file_name][2]
						tissue_type = dictionary_total_reads_bam_files[bam_file_name][3]
						short_list = dictionary_total_reads_bam_files[bam_file_name][4]

						if tissue_name == '':
							print (bam_file_name, 'There is a problem to define the tissue for this bamfile when analyzing this peptide : ', peptide)
						
						key = (tissue_name, tissue_type, short_list)
						try:
							info_tissue = info_tissue_peptide[key]
							try:
								info_tissue[peptide].append(norm)
							except KeyError:
								info_tissue[peptide] = [norm]
						except KeyError:
							info_tissue_peptide[key] = {peptide: [norm]}

				for key_info_tissue, peptides in info_tissue_peptide.items():
					
					tissue_name = key_info_tissue[0]
					tissue_type = key_info_tissue[1]
					short_list = key_info_tissue[2]

					for peptide, info in peptides.items():
						mean = np.mean(info)
						median = np.median(info)
						# This is to get for each peptide the normalisation even though the peptide is several peptide types.
						# Independent normalisation for the same peptide in each peptide type
						# peptide_type = self.peptides_types[peptide].split(';')
						# for type_ in peptide_type:
						# 	aux = [peptide, type_, tissue_name, tissue_type, short_list, median, mean]
						# 	data.append(aux)
						peptide_type = self.peptides_types[peptide]
						aux = [peptide, peptide_type, tissue_name, tissue_type, short_list, median, mean]
						data.append(aux)

				df_norm = pd.DataFrame(data, columns=['Peptide','Peptide_Type', 'Tissue', 'Tissue_type', 'Short_list','median','mean'])
				df_norm.to_csv(self.path_to_output_total_transcription_expression_heatmap_folder+'norm_info.csv', index=False)

			def_norm_log.reset_index(inplace=True)
			
			def_norm_log.to_csv(self.path_to_output_temps_folder+self.name_exp+type_save, index=False, header=True)
			self.super_logger.info('Normalization Information saved to : %s ', self.path_to_output_temps_folder+self.name_exp+type_save)

			t_2 = time.time()
			total = t_2-t_0
			self.super_logger.info('Total time run function get_normalization to end : %f min', (total/60.0))
		else:
			self.super_logger.info('Normalization information already collected in the output folder : %s --> Skipping this step!', self.path_to_output_temps_folder+self.name_exp+type_save)
			def_norm_log = pd.read_csv(self.path_to_output_temps_folder+self.name_exp+type_save)
		
		return def_norm_log

		
# 				t1 				ts2 			ts3 			ts4
# rphm			x				x				x				x
# log(rphm)		log(x+1,10)		log(x+1,10)		log(x+1,10)		log(x+1,10)
# mean_log		sum(log(x+1,10)	+ log(x+1,10)+log(x+1,10)+log(x+1,10))/4
# arithmetic_mean	sum(x1	+ x2 + x3 + x4)/4

