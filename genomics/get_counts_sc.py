import warnings, os
warnings.filterwarnings("ignore")
import os, time, subprocess, pickle, multiprocessing, os, _thread, csv, collections, pysam, copy
import genomics.get_alignments as get_alig
import pandas as pd
from pathos.multiprocessing import ProcessPool
import utils.useful_functions as uf

__author__ = "Maria Virginia Ruiz Cuevas"
__email__ = "maria.virginia.ruiz.cuevas@umontreal.ca"


NUM_WORKERS =  multiprocessing.cpu_count()

path_to_lib = '/'.join(os.path.abspath(__file__).split('/')[:-3])+'/lib/'

class GetCountsSC:

	def __init__(self, path_to_output_folder, name_exp, mode, light, peptides_by_type, super_logger):
		self.light = light
		if self.light:
			self.path_to_output_folder = path_to_output_folder+'res_light/'
		else:
			self.path_to_output_folder = path_to_output_folder+'res/'
		
		self.name_exp = name_exp
		self.mode = mode
		self.path_to_output_folder_alignments = path_to_output_folder+'alignments/'
		self.peptides_by_type = peptides_by_type
		self.super_logger = super_logger
		

	def get_counts(self, perfect_alignments, bam_files_list):

		df_counts = pd.DataFrame()
		df_counts_filtered = pd.DataFrame()

		exist_rna = os.path.exists(self.path_to_output_folder+self.name_exp+'_rna_count.csv')
		
		last_treated_bam_file = os.path.exists(self.path_to_output_folder_alignments+'info_trated_bam_files.pkl')

		if last_treated_bam_file:

			with open(self.path_to_output_folder_alignments+'info_trated_bam_files.pkl', 'rb') as fp:
				last_treated_bam_file = pickle.load(fp)

			try:
				with open(self.path_to_output_folder_alignments+'peptides_info.pkl', 'rb') as fp:
					peptides_info = pickle.load(fp)

			except FileNotFoundError:
				pass

		else:
			last_treated_bam_file = -1
		
		
		if self.light:
			alignment_information = self.path_to_output_folder_alignments+'Alignments_information_light.dic'
		else:
			alignment_information = self.path_to_output_folder_alignments+'Alignments_information.dic'
			
		
		if not exist_rna:
			t_0 = time.time()
			
			info_bams = []
			bams = []
			
			list_bam_files_order = []
			index_sample = 0
			for name_sample, info_bam in sorted(bam_files_list.items(), key=lambda e: e[1][-1], reverse=False):
				info_bams.append((index_sample,info_bam))
				list_bam_files_order.append(name_sample)
				index_sample += 1
			
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

			peptides_info = {}
			cell_lines = set()

			for idx, bam_file in enumerate(info_bams):
				
				if idx > last_treated_bam_file:

					t0_bam_file = time.time()
					bams = [bam_file] * len(keys)
					results = pool.map(self.get_counts_sample, bams, keys, values)
					not_permission = False

					for res in results:
						for index, count_align in enumerate(res):

							if index == 0: 
								peptide = count_align[0]
								alignment = count_align[1]
								index_sample = count_align[2]
								strand = count_align[3]
								key = alignment+'_'+strand
								try:
									peptide_info_aux = peptides_info[peptide]
									try:
										peptide_info_aux_alignment = peptide_info_aux[key]
									except KeyError:
										peptide_info_aux[key] = {}
								except KeyError:
									peptides_info[peptide] = {key: {}}
								
							else:
								count_info = count_align[0]
								if count_info == -1:
									not_permission = True
								sequence = count_align[1]
								
								if len(count_info) > 0:
									cells = set(count_info.keys())
									cell_lines = cell_lines.union(cells)

								info_alignment = peptides_info[peptide][key]
								try:
									info_sequence = info_alignment[sequence]
									if len(count_info) > 0:
										info_sequence.update(count_info)

								except KeyError:
									info_alignment[sequence] = count_info
									
								count = sum(count_info.values())
								new_key = peptide+'_'+alignment+'_'+sequence
								if len(perfect_alignments[new_key][-2]) == 0:
									perfect_alignments[new_key][-2] = [0]*total_samples
									
								perfect_alignments[new_key][-2][index_sample] = count


					t1_bam_file = time.time()
					time_final = (t1_bam_file-t0_bam_file)/60.0

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

						with open(self.path_to_output_folder_alignments+'peptides_info.pkl', 'wb') as handle:
							pickle.dump(peptides_info, handle, protocol=pickle.HIGHEST_PROTOCOL)
				else:
					self.super_logger.info('Bam File already processed: %s %s.', str(idx), bam_file[0])
			
			pool.close()
			pool.join()
			pool.clear()

			with open(self.path_to_output_folder_alignments+'info_trated_bam_files.pkl', 'wb') as f:  
				pickle.dump(idx, f)
			
			with open(alignment_information, 'wb') as handle:
				pickle.dump(perfect_alignments, handle, protocol=pickle.HIGHEST_PROTOCOL)

			peptides_order = []
			
			for peptide_type, peptides in self.peptides_by_type.items():
				peptides_order.extend(peptides)

			cell_lines = list(cell_lines)

			header = ['Peptide', 'Alignment', 'MCS', 'Strand']+cell_lines
			data = []

			for peptide in peptides_order:

				try:
					info_peptide = peptides_info[peptide]
					
					for key, info_sequences in info_peptide.items():
						key_splited = key.split('_')
						alignment = key_splited[0]
						strand = key_splited[1]
						
						for sequence, value_sequence in info_sequences.items():
							if len(value_sequence)>0:
								aux = [peptide, alignment, sequence, strand]
								zeros = [0]*len(cell_lines)
								for cell, count in value_sequence.items():
									index = cell_lines.index(cell)
									zeros[index] = count
								aux = aux+zeros
								data.append(aux)

				except KeyError:
					pass


			df_alignments = pd.DataFrame(data, columns = header)
			df_counts = df_alignments.groupby(['Peptide']).sum().reset_index()
			
			sorterIndex = dict(zip(peptides_order, range(len(peptides_order))))
			df_counts['order'] = df_counts['Peptide'].map(sorterIndex)
			df_counts.sort_values(['order'], ascending = [True], inplace = True)
			df_counts.drop('order', 1, inplace = True)
			df_counts.set_index('Peptide',inplace=True)
			
			name_path = self.path_to_output_folder+self.name_exp+'_rna_count_All_alignments.csv'
			df_alignments.to_csv(name_path, index=False, header=True)

			df_counts.to_csv(self.path_to_output_folder+self.name_exp+'_rna_count.csv', index=True, header=True)
			self.super_logger.info('Counts Information saved to : %s ', self.path_to_output_folder+self.name_exp+'_rna_count.csv')

			t_2 = time.time()
			total = t_2-t_0
			self.super_logger.info('Total time run function get_counts to end : %f min', (total/60.0))
		else:
			self.super_logger.info('Count information already collected in the output folder : %s --> Skipping this step!', self.path_to_output_folder+self.name_exp+'_rna_count.csv')
			df_counts = pd.read_csv(self.path_to_output_folder+self.name_exp+'_rna_count.csv', index_col=0)
			
			with open(alignment_information, 'rb') as fp:
				perfect_alignments = pickle.load(fp)

			df_alignments = pd.read_csv(self.path_to_output_folder+self.name_exp+'_rna_count_All_alignments.csv', index_col=0)

		return df_counts, perfect_alignments, df_alignments


	def get_counts_sample(self, bam, peptide_alignment, sequences):
		
		index_sample = bam[0]
		bam_file = bam[1][0]
		library = bam[1][1]
		sens = bam[1][2]

		count = 0
		peptide = peptide_alignment.split('_')[0]
		alignment = peptide_alignment.split('_')[1]
		strand = peptide_alignment.split('_')[2]

		chr = alignment.split(':')[0]
		
		region_to_query = chr+':'+alignment.split(':')[1].split('-')[0]+'-'+alignment.split(':')[1].split('-')[-1]

		contReads_to_return = self.get_depth_with_view(region_to_query, bam_file, index_sample, library, sens, strand, sequences)

		to_return = [[peptide, alignment, index_sample, strand]]
		
		for index, sequence in enumerate(sequences):
			try:
				count = contReads_to_return[sequence]
			except IndexError:
				count = -1
			to_return.append([count, sequence])
		 
		return to_return


	def get_depth_with_view(self, region_to_query, bam_file, index_sample, library, sens, strand, sequences):

		contReads_to_return = {}

		library = library.lower()
		sens = sens.lower()
		
		if library == 'unstranded':
			try:
				count_1 = pysam.view("-F0X100", bam_file, region_to_query)
			except pysam.utils.SamtoolsError: 
				return -1

			list_seq = set()
			for seq in sequences:
				contReads = 0
				rcmcs = uf.reverseComplement(seq)
				contReads += count_1.count(seq)
				contReads += count_1.count(rcmcs)
				contReads_to_return[seq] = {}
				if contReads > 0:
					list_seq.add(seq)

			if len(list_seq) > 0:
				for seq in list_seq:
					rcmcs = uf.reverseComplement(seq)
					for read in count_1.split('\n'):
						if seq in read :
							try:
								cell = read.split('CB:Z:')[1].split('-')[0]
								name_cell = str(index_sample)+'_'+cell
								try:
									contReads_to_return[seq][name_cell] += 1
								except KeyError:
									contReads_to_return[seq][name_cell] = 1
							except IndexError:
								try:
									cell = read.split('CR:Z:')[1].split('\t')[0]
									name_cell = str(index_sample)+'_'+cell
									try:
										contReads_to_return[seq][name_cell] += 1
									except KeyError:
										contReads_to_return[seq][name_cell] = 1
								except IndexError:
									pass

						if rcmcs in read :
							try:
								cell = read.split('CB:Z:')[1].split('-')[0]
								name_cell = str(index_sample)+'_'+cell
								try:
									contReads_to_return[seq][name_cell] += 1
								except KeyError:
									contReads_to_return[seq][name_cell] = 1
							except IndexError:
								try:
									cell = read.split('CR:Z:')[1].split('\t')[0]
									name_cell = str(index_sample)+'_'+cell
									try:
										contReads_to_return[seq][name_cell] += 1
									except KeyError:
										contReads_to_return[seq][name_cell] = 1
								except IndexError:
									pass

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

			list_seq = set()
			for seq in sequences:
				contReads = count_1.count(seq)
				contReads_to_return[seq] = {}
				if contReads > 0:
					list_seq.add(seq)

			if len(list_seq) > 0:
				for seq in list_seq:
					for read in count_1.split('\n'):
						if seq in read:
							try:
								cell = read.split('CB:Z:')[1].split('-')[0]
								name_cell = str(index_sample)+'_'+cell
								try:
									contReads_to_return[seq][name_cell] += 1
								except KeyError:
									contReads_to_return[seq][name_cell] = 1
							except IndexError:
								try:
									cell = read.split('CR:Z:')[1].split('\t')[0]
									name_cell = str(index_sample)+'_'+cell
									try:
										contReads_to_return[seq][name_cell] += 1
									except KeyError:
										contReads_to_return[seq][name_cell] = 1
								except IndexError:
									pass
							
				
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
				contReads_to_return[seq] = {}

				for read in count_1_split:
					if seq in read :
						name = read.split('\t')[0]
						try:
							cell = read.split('CB:Z:')[1].split('-')[0]
							name_cell = str(index_sample)+'_'+cell
							if name not in reads_name:
								reads_name.add(name)
								contReads_to_return[seq][name_cell] += 1
						except IndexError:
							try:
								cell = read.split('CR:Z:')[1].split('\t')[0]
								name_cell = str(index_sample)+'_'+cell
								if name not in reads_name:
									reads_name.add(name)
									contReads_to_return[seq][name_cell] += 1
							except IndexError:
								pass

				for read in count_2_split:
					if rcmcs in read :
						name = read.split('\t')[0]
						try:
							cell = read.split('CB:Z:')[1].split('-')[0]
							name_cell = str(index_sample)+'_'+cell
							if name not in reads_name:
								reads_name.add(name)
								contReads_to_return[seq][name_cell] += 1
						except IndexError:
							try:
								cell = read.split('CR:Z:')[1].split('\t')[0]
								name_cell = str(index_sample)+'_'+cell
								if name not in reads_name:
									reads_name.add(name)
									contReads_to_return[seq][name_cell] += 1
							except IndexError:
								pass

			# reading htslib https://github.com/DecodeGenetics/graphtyper/issues/57
			# From https://davetang.org/muse/2018/06/06/10x-single-cell-bam-files/ I know the CR and CB differences
			# https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/bam
		return contReads_to_return

	
