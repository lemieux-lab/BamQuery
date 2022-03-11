import warnings, os
warnings.filterwarnings("ignore")
import os, time, subprocess, pickle, multiprocessing, os, _thread, csv, collections, pysam, copy
import genomics.get_counts_from_sample as get_counts_sample
from readers.intersection_alignments_annotations import IntersectAnnotations
from genomics.get_information_from_bed_intersection import GetInformationBEDIntersection
import pandas as pd
from pathos.multiprocessing import ProcessPool
import utils.useful_functions as uf
import numpy as np
import math

__author__ = "Maria Virginia Ruiz Cuevas"
__email__ = "maria.virginia.ruiz.cuevas@umontreal.ca"


NUM_WORKERS =  multiprocessing.cpu_count()

path_to_lib = '/'.join(os.path.abspath(__file__).split('/')[:-3])+'/lib/'

class GetCounts:

	def __init__(self, path_to_output_folder, name_exp, mode, light, peptides_by_type, super_logger):
		self.light = light
		self.mode = mode
		self.path_to_output_folder = path_to_output_folder
		self.path_to_output_folder_alignments = path_to_output_folder+'alignments/'
		self.name_exp = name_exp
		self.peptides_by_type = peptides_by_type
		self.super_logger = super_logger
		
		if self.light:
			self.path_to_output_temps_folder = path_to_output_folder+'res_light/temps_files/'
			self.alignment_information_path = self.path_to_output_folder_alignments+'Alignments_information_light_rna.dic'

		else:
			self.path_to_output_temps_folder = path_to_output_folder+'res/temps_files/'
			self.alignment_information_path = self.path_to_output_folder_alignments+'Alignments_information_rna.dic'
		
		self.count_path = self.path_to_output_temps_folder+self.name_exp+'_rna_count.csv'
		self.count_all_alignments_path = self.path_to_output_temps_folder+self.name_exp+'_rna_count_All_alignments.csv'
		self.last_treated_bam_file = self.path_to_output_folder_alignments+'info_treated_bam_files.pkl'
		self.to_write = self.path_to_output_folder_alignments+'to_write.pkl'

		
		if self.mode == 'translation':
			self.path_to_output_temps_folder = path_to_output_folder+'res_translation/temps_files/'
			self.path_to_output_bed_folder = path_to_output_folder+'res_translation/BED_files/'
			self.count_path = self.path_to_output_temps_folder+self.name_exp+'_ribo_counts.csv'
			self.count_all_alignments_path = self.path_to_output_temps_folder+self.name_exp+'_ribo_count_All_alignments.csv'
			self.last_treated_bam_file = self.path_to_output_folder_alignments+'info_treated_ribo_bam_files.pkl'
			self.to_write = self.path_to_output_folder_alignments+'to_ribo_write.pkl'
			if self.light:
				self.alignment_information_path = self.path_to_output_folder_alignments+'Alignments_information_light_ribo.dic'
			else:
				self.alignment_information_path = self.path_to_output_folder_alignments+'Alignments_information_ribo.dic'

		

	def get_coverage(self, perfect_alignments, bam_files_list, genome_version):

		times = []
		df_counts = pd.DataFrame()
			
		ribo_coverage_info = self.path_to_output_temps_folder+self.name_exp+'_ribo_coverage_info.csv'
		ribo_coverage_info_all_aligments_path = self.path_to_output_temps_folder+self.name_exp+'_ribo_coverage_info_All_alignments.csv'
		exists_ribo_coverage_info = os.path.exists(ribo_coverage_info)
		
		alignment_information_ribo_path = self.path_to_output_folder_alignments+'Alignments_information_ribo_assembly.dic'
		exists_alignment_information_ribo = os.path.exists(alignment_information_ribo_path)
		
		last_treated_assembly = os.path.exists(self.path_to_output_folder_alignments+'info_treated_assembly.pkl')
		if last_treated_assembly:

			with open(self.path_to_output_folder_alignments+'info_treated_assembly.pkl', 'rb') as fp:
				last_treated_assembly = pickle.load(fp)

			try:
				with open(self.path_to_output_folder_alignments+'to_write.pkl', 'rb') as fp:
					to_write = pickle.load(fp)
			except FileNotFoundError:
				to_write = {} 
		else:
			last_treated_assembly = -1
			to_write = {} 


		if  not exists_ribo_coverage_info :
			t_0 = time.time()
			
			if not exists_alignment_information_ribo :
				alignment_information_ribo = self.transform_initial_alignment_information(perfect_alignments, alignment_information_ribo_path)
			else:
				with open(alignment_information_ribo_path, 'rb') as fp:
			 		alignment_information_ribo = pickle.load(fp)

			list_bam_files_order = []
			for name_sample, info_bam in sorted(bam_files_list.items(), key=lambda e: e[1][-2], reverse=False):
				list_bam_files_order.append(name_sample)
			
			def get_index(sample):
				index = list_bam_files_order.index(sample)
				return index + 1

			total_samples = len(list_bam_files_order)

			intersect_to_annotations = IntersectAnnotations(alignment_information_ribo, self.path_to_output_folder, self.mode, self.name_exp, self.super_logger, genome_version)
			exists_bed_to_intersect = os.path.exists(self.path_to_output_bed_folder+'to_intersect_to_annotations.bed')
			
			if not exists_bed_to_intersect:
				intersect_to_annotations.generate_BED_files()

			files_to_intersect = []

			for bam_file_name, bam_file_info in bam_files_list.items():
				bed_intersected = bam_file_info[2]
				intersected_bed = os.path.exists(bed_intersected)
				if not intersected_bed :
					files_to_intersect.append((bam_file_info[1], bed_intersected))

			if len(files_to_intersect) > 0:
				intersect_to_annotations.perform_intersection_with_assemblies(files_to_intersect)

			information_from_bed_files_intersected = self.path_to_output_bed_folder+'/information_from_bed_files_intersected.dic'
			information_transcripts_intersected = self.path_to_output_bed_folder+'/transcripts_intersected.dic'

			resume_bed_file_intersected = os.path.exists(information_from_bed_files_intersected)
			resume_transcripts = os.path.exists(information_transcripts_intersected)
			
			peptides_translated = {}
			self.transcripts_intersected = {}

			if not resume_bed_file_intersected or not resume_transcripts:
				intersected_bed_files_translation = GetInformationBEDIntersection(self.path_to_output_folder, self.mode)
				intersected_bed_files_translation.get_ribosome_profiling_transcripts_overlap(bam_files_list)
				peptides_translated = intersected_bed_files_translation.peptides_translated
				self.transcripts_intersected = intersected_bed_files_translation.transcripts_intersected
			else:
				with open(information_from_bed_files_intersected, 'rb') as fp:
					peptides_translated = pickle.load(fp)
				with open(information_transcripts_intersected, 'rb') as fp:
					self.transcripts_intersected = pickle.load(fp)
			
			pool = ProcessPool(nodes=NUM_WORKERS)

			for idx, bam_file in enumerate(bam_files_list):
				
				if idx > last_treated_assembly:

					t0_bam_file = time.time()

					keys = list(peptides_translated[bam_file].keys())
					values = list(peptides_translated[bam_file].values())
					
					results = pool.map(self.get_tpm_for_overlap_peptide, keys, values, [bam_file]*len(keys))	
					
					for res in results:
						key = res[0]
						peptide = key.split('_')[0]
						tpm_mean = res[1]

						if len(alignment_information_ribo[key]) == 1:
							alignment_information_ribo[key].append([0]*total_samples)
						
						index = get_index(bam_file) 
						alignment_information_ribo[key][1][index-1] = tpm_mean
						strand = alignment_information_ribo[key][0]

						try:
							peptide_info_to_write = to_write[peptide]

							try:
								peptide_info_to_write[key][index] = tpm_mean
							except KeyError:
								counts = [0]*total_samples
								to_add = [strand]
								to_add.extend(counts)
								to_add[index] = tpm_mean
								peptide_info_to_write[key] = to_add

						except KeyError:
							counts = [0]*total_samples
							to_add = [strand]
							to_add.extend(counts)
							to_add[index]  = tpm_mean
							to_write[peptide] = {key: to_add}


					t1_bam_file = time.time()
					time_final = (t1_bam_file-t0_bam_file)/60.0
					times.append(time_final)

					if (idx % 100 == 0) and (idx != 0):
						self.super_logger.info('Saving information for Bam Files processed')

						with open(self.path_to_output_folder_alignments+'info_treated_assembly.pkl', 'wb') as f:  
							pickle.dump(idx, f)

						with open(alignment_information_ribo_path, 'wb') as handle:
							pickle.dump(alignment_information_ribo, handle, protocol=pickle.HIGHEST_PROTOCOL)

						with open(self.path_to_output_folder_alignments+'to_write.pkl', 'wb') as handle:
							pickle.dump(to_write, handle, protocol=pickle.HIGHEST_PROTOCOL)
			
				else:
					self.super_logger.info('Bam File already processed: %s %s.', str(idx), bam_file)

			pool.close()
			pool.join()
			pool.clear()

			self.super_logger.info('Average time to process a BamFile : %f min', np.mean(times))

			with open(self.path_to_output_folder_alignments+'info_treated_assembly.pkl', 'wb') as f:  
				pickle.dump(idx, f)
			
			with open(alignment_information_ribo_path, 'wb') as handle:
				pickle.dump(alignment_information_ribo, handle, protocol=pickle.HIGHEST_PROTOCOL)

			header = ['Peptide Type','Peptide', 'Position', 'Strand']
			header.extend(list_bam_files_order)
			to_write_list = []
			
			for peptide, info_peptide in to_write.items():
				to_add = []
				peptide_type = self.peptides_by_type[peptide]
				for alignment, info_counts in info_peptide.items():
					position = alignment.split('_')[1]
					to_add =  [peptide_type, peptide, position]
					to_add.extend(info_counts)
					to_write_list.append(to_add)

			df_alignments = pd.DataFrame(to_write_list, columns=header)
			df_counts = df_alignments.groupby(['Peptide Type', 'Peptide']).mean().reset_index()
			df_counts.sort_values(by=['Peptide Type'])
			df_counts.set_index(['Peptide Type', 'Peptide'], inplace=True)
			df_counts.to_csv(ribo_coverage_info,  header=True)
			df_alignments.to_csv(ribo_coverage_info_all_aligments_path, index=False, header=True)

			self.super_logger.info('Coverage Information saved to : %s ', ribo_coverage_info)

			t_2 = time.time()
			total = t_2-t_0
			self.super_logger.info('Total time run function get_counts to end : %f min', (total/60.0))

		else:
			self.super_logger.info('TPM information already collected in the output folder : %s --> Skipping this step!', ribo_coverage_info)
			df_counts = pd.read_csv(ribo_coverage_info, index_col=[0,1])
			
			with open(alignment_information_ribo_path, 'rb') as fp:
				alignment_information_ribo = pickle.load(fp)

			df_alignments = pd.read_csv(ribo_coverage_info_all_aligments_path, index_col=[0,1])

		
		return df_counts, alignment_information_ribo, df_alignments

	
	def transform_initial_alignment_information(self, perfect_alignments, path_alignment_information_ribo):

		alignment_information_ribo = {}
		for peptide_key, information_peptide in perfect_alignments.items():
			peptide = peptide_key.split('_')[0]
			if peptide in self.peptides_by_type :
				new_key =  '_'.join(peptide_key.split('_')[:-1])
				alignment_information_ribo[new_key] = [information_peptide[0]]

		with open(path_alignment_information_ribo, 'wb') as handle:
			pickle.dump(alignment_information_ribo, handle, protocol=pickle.HIGHEST_PROTOCOL)

		return alignment_information_ribo


	def get_tpm_for_overlap_peptide(self, key, values, sample):

		def get_tpm(transcript):
			return self.transcripts_intersected[transcript+'_'+sample]+1

		values = list(values.values())
		transcripts_intersected = values[0].intersection(*values)
		tpm_mean = np.mean(np.log10(list(map(get_tpm, transcripts_intersected))))
		return key, tpm_mean


	def get_counts(self, perfect_alignments, bam_files_list, overlap):

		times = []
		df_counts = pd.DataFrame()

		exists_count = os.path.exists(self.count_path)
		exists_alignment_information = os.path.exists(self.alignment_information_path)
		
		last_treated_bam_file = os.path.exists(self.last_treated_bam_file )
		if last_treated_bam_file:

			with open(self.last_treated_bam_file , 'rb') as fp:
				last_treated_bam_file = pickle.load(fp)

			try:
				with open(self.to_write, 'rb') as fp:
					to_write = pickle.load(fp)
			except FileNotFoundError:
				to_write = {} 
		else:
			last_treated_bam_file = -1
			to_write = {} 
		

		if not exists_count :
			t_0 = time.time()
			
			if not exists_alignment_information :
				alignment_information = copy.deepcopy(perfect_alignments)
			else:
				with open(self.alignment_information_path, 'rb') as fp:
			 		alignment_information = pickle.load(fp)

			info_bams = []
			bams = []

			list_bam_files_order = []
			for name_sample, info_bam in sorted(bam_files_list.items(), key=lambda e: e[1][-2], reverse=False):
				info_bams.append((name_sample,info_bam))
				list_bam_files_order.append(name_sample)
			
			def get_index(sample):
				index = list_bam_files_order.index(sample)
				return index + 1

			total_samples = len(list_bam_files_order)

			keys = alignment_information.keys()
			self.super_logger.info('Total MCS mapped : %s ', str(len(keys)))

			modif_dic = {}
			
			for key in keys:
				split_key = key.split('_')
				peptide = split_key[0]
				position = split_key[1]
				seq = split_key[2]

				strand = alignment_information[key][0]
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
					results = pool.map(get_counts_sample.get_counts_sample, bams, keys, values, [overlap]*len(values))
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
								if len(alignment_information[key][-1]) == 0:
									alignment_information[key][-1] = [0]*total_samples

								alignment_information[key][-1][index-1] = count


					t1_bam_file = time.time()
					time_final = (t1_bam_file-t0_bam_file)/60.0
					times.append(time_final)

					if not_permission:
						self.super_logger.info('Bam File : %s %s couldn\'t be processed. Failed to open, permission denied. Time : %f min', str(idx), bam_file[0], time_final)
					else:
						self.super_logger.info('Processed Bam File : %s %s. Time : %f min', str(idx), bam_file[0], time_final)

					if (idx % 100 == 0) and (idx != 0):
						self.super_logger.info('Saving information for Bam Files processed')

						with open(self.last_treated_bam_file, 'wb') as f:  
							pickle.dump(idx, f)

						with open(self.alignment_information_path, 'wb') as handle:
							pickle.dump(alignment_information, handle, protocol=pickle.HIGHEST_PROTOCOL)

						with open(self.to_write, 'wb') as handle:
							pickle.dump(to_write, handle, protocol=pickle.HIGHEST_PROTOCOL)
						
				else:
					self.super_logger.info('Bam File already processed: %s %s.', str(idx), bam_file[0])

			pool.close()
			pool.join()
			pool.clear()

			self.super_logger.info('Average time to process a BamFile : %f min', np.mean(times))

			with open(self.last_treated_bam_file, 'wb') as f:  
				pickle.dump(idx, f)
			
			with open(self.alignment_information_path, 'wb') as handle:
				pickle.dump(alignment_information, handle, protocol=pickle.HIGHEST_PROTOCOL)

			header = ['Peptide Type', 'Peptide', 'Position', 'MCS', 'Strand']
			header.extend(list_bam_files_order)
			to_write_list = []
			
			for peptide, info_peptide in to_write.items():
				to_add = []
				peptide_type = self.peptides_by_type[peptide]
				for alignment, info_counts in info_peptide.items():
					position = alignment.split('_')[0]
					MCS = alignment.split('_')[1]
					to_add =  [peptide_type, peptide, position, MCS]
					to_add.extend(info_counts)
					to_write_list.append(to_add)

			df_alignments = pd.DataFrame(to_write_list, columns=header)
			df_counts = df_alignments.groupby(['Peptide Type', 'Peptide']).sum().reset_index()
			df_counts.sort_values(by=['Peptide Type'])
			df_counts.set_index(['Peptide Type', 'Peptide'], inplace=True)
			df_counts.to_csv(self.count_path, header=True)
			
			df_alignments.to_csv(self.count_all_alignments_path, index=False, header=True)

			self.super_logger.info('Counts Information saved to : %s ', self.count_path)

			t_2 = time.time()
			total = t_2-t_0
			self.super_logger.info('Total time run function get_counts to end : %f min', (total/60.0))

		else:
			self.super_logger.info('Count information already collected in the output folder : %s --> Skipping this step!', rna_count_path)
			df_counts = pd.read_csv(self.count_path, index_col=[0,1])

			with open(self.alignment_information_path, 'rb') as fp:
				alignment_information = pickle.load(fp)

			df_alignments = pd.read_csv(self.count_all_alignments_path, index_col=[0,1])

		try:
			os.remove(self.to_write)
		except FileNotFoundError:
			pass
		
		return df_counts, alignment_information, df_alignments




	
