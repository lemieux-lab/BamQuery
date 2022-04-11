import warnings, os
warnings.filterwarnings("ignore")
import os, time, subprocess, pickle, multiprocessing, os, _thread, csv, collections, pysam, copy
import genomics.get_alignments as get_alig
import pandas as pd
from pathos.multiprocessing import ProcessPool
import utils.useful_functions as uf
from itertools import repeat
from operator import itemgetter
import re

__author__ = "Maria Virginia Ruiz Cuevas"
__email__ = "maria.virginia.ruiz.cuevas@umontreal.ca"


NUM_WORKERS =  multiprocessing.cpu_count()

path_to_lib = '/'.join(os.path.abspath(__file__).split('/')[:-3])+'/lib/'

class GetCountsSC:

	def __init__(self, path_to_output_folder, name_exp, mode, light, peptides_by_type, super_logger):
		self.path_to_output_folder = path_to_output_folder+'res/'
		self.name_exp = name_exp
		self.mode = mode
		self.path_to_output_folder_alignments = path_to_output_folder+'alignments/'
		self.peptides_by_type = peptides_by_type
		self.super_logger = super_logger
	

	def filter_alignment_information(self, perfect_alignments, path_alignment_information):

		alignment_information = {}
		for peptide_key, information_peptide in perfect_alignments.items():
			peptide = peptide_key.split('_')[0]
			if peptide in self.peptides_by_type :
				alignment_information[peptide_key] = information_peptide

		with open(path_alignment_information, 'wb') as handle:
			pickle.dump(alignment_information, handle, protocol=pickle.HIGHEST_PROTOCOL)

		return alignment_information


	def get_counts(self, perfect_alignments, bam_files_list):


		def sum_dict(d1, d2):
			for key, value in d1.items():
				d1[key] = value + d2.get(key, 0)
			return d1
	
		df_counts = pd.DataFrame()
		
		rna_sc_count_path = self.path_to_output_folder+self.name_exp+'_rna_sc_count.csv'
		rna_sc_count_all_alignments_path = self.path_to_output_folder+self.name_exp+'_rna_sc_count_All_alignments.csv'
		
		exists_rna_sc = os.path.exists(rna_sc_count_path)
		
		alignment_information_sc_path = self.path_to_output_folder_alignments+'Alignments_information_sc.dic'
		exists_alignment_information_sc = os.path.exists(alignment_information_sc_path)
		
		last_treated_bam_file = os.path.exists(self.path_to_output_folder_alignments+'info_treated_bam_files.pkl')
		
		if last_treated_bam_file:

			with open(self.path_to_output_folder_alignments+'info_trated_bam_files.pkl', 'rb') as fp:
				last_treated_bam_file = pickle.load(fp)

			try:
				with open(self.path_to_output_folder_alignments+'peptides_info.pkl', 'rb') as fp:
					peptides_info = pickle.load(fp)

			except FileNotFoundError:
				peptides_info = {}

		else:
			last_treated_bam_file = -1
			peptides_info = {}
		

		if not exists_rna_sc:
			t_0 = time.time()
			
			if not exists_alignment_information_sc :
				alignment_information_sc = self.filter_alignment_information(perfect_alignments, alignment_information_sc_path)
			else:
				with open(alignment_information_sc_path, 'rb') as fp:
			 		alignment_information_sc = pickle.load(fp)

			info_bams = []
			bams = []
			
			list_bam_files_order = []
			
			index_sample = 0
			for name_sample, info_bam in sorted(bam_files_list.items(), key=lambda e: e[1][-2], reverse=False):
				info_bams.append((index_sample,info_bam))
				list_bam_files_order.append(name_sample)
				index_sample += 1

			def get_index(sample):
				index = list_bam_files_order.index(sample)
				return index + 1

			total_samples = len(list_bam_files_order)

			keys = alignment_information_sc.keys()
			self.super_logger.info('Total MCS mapped : %s ', str(len(keys)))

			modif_dic = {}
			for key in keys:
				split_key = key.split('_')
				peptide = split_key[0]
				position = split_key[1]
				seq = split_key[2]

				strand = alignment_information_sc[key][0]
				new_key = peptide+'_'+position+'_'+strand
				
				try:
					modif_dic[new_key].append(seq)
				except KeyError:
					modif_dic[new_key] = [seq]

			keys = list(modif_dic.keys())
			values = list(modif_dic.values())
			self.super_logger.info('Total unique regions : %s ', str(len(keys)))

			pool = ProcessPool(nodes=NUM_WORKERS)

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
										dict1 = [info_sequence, count_info]
										info_alignment[sequence] = reduce(sum_dict, dict1)
								
								except KeyError:
									info_alignment[sequence] = count_info

								count = sum(count_info.values())
								new_key = peptide+'_'+alignment+'_'+sequence
								if len(alignment_information_sc[new_key][-1]) == 0:
									alignment_information_sc[new_key][-1] = [0]*total_samples
									
								alignment_information_sc[new_key][-1][index_sample] = count


					t1_bam_file = time.time()
					time_final = (t1_bam_file-t0_bam_file)/60.0

					if not_permission:
						self.super_logger.info('Bam File : %s %s couldn\'t be processed. Failed to open, permission denied. Time : %f min', str(idx), bam_file[0], time_final)
					else:
						self.super_logger.info('Processed Bam File : %s %s. Time : %f min', str(idx), list_bam_files_order[idx], time_final)

					if (idx % 100 == 0) and (idx != 0):
						self.super_logger.info('Saving information for Bam Files processed')

						with open(self.path_to_output_folder_alignments+'info_treated_bam_files.pkl', 'wb') as f:  
							pickle.dump(idx, f)

						with open(alignment_information_sc_path, 'wb') as handle:
							pickle.dump(alignment_information_sc, handle, protocol=pickle.HIGHEST_PROTOCOL)

						with open(self.path_to_output_folder_alignments+'peptides_info.pkl', 'wb') as handle:
							pickle.dump(peptides_info, handle, protocol=pickle.HIGHEST_PROTOCOL)
				else:
					self.super_logger.info('Bam File already processed: %s %s.', str(idx), bam_file[0])
			
			pool.close()
			pool.join()
			pool.clear()

			with open(self.path_to_output_folder_alignments+'info_trated_bam_files.pkl', 'wb') as f:  
				pickle.dump(idx, f)
			
			with open(alignment_information_sc_path, 'wb') as handle:
				pickle.dump(alignment_information_sc, handle, protocol=pickle.HIGHEST_PROTOCOL)

			cell_lines = list(cell_lines)
			header = ['Peptide Type', 'Peptide', 'Position', 'MCS', 'Strand']+cell_lines
			data = []
			
			for peptide, info_peptide in peptides_info.items():
				to_add = []
				peptide_type = self.peptides_by_type[peptide]
				for alignment, info_sequences in info_peptide.items():
					key_splited = alignment.split('_')
					alignment = key_splited[0]
					strand = key_splited[1]
					
					for sequence, value_sequence in info_sequences.items():
						if len(value_sequence)>0:
							aux = [peptide_type, peptide, alignment, sequence, strand]
							zeros = [0]*len(cell_lines)
							for cell, count in value_sequence.items():
								index = cell_lines.index(cell)
								zeros[index] = count
							aux = aux+zeros
							data.append(aux)

			df_alignments = pd.DataFrame(data, columns = header)
			df_counts = df_alignments.groupby(['Peptide Type', 'Peptide']).sum().reset_index()
			df_counts.sort_values(by=['Peptide Type'])

			df_counts.to_csv(rna_sc_count_path,  header=True, index=False)
			df_alignments.to_csv(rna_sc_count_all_alignments_path, index=False, header=True)

			self.super_logger.info('Counts SC Information saved to : %s ', rna_sc_count_path)

			t_2 = time.time()
			total = t_2-t_0
			self.super_logger.info('Total time run function get_counts to end : %f min', (total/60.0))
			
		else:
			self.super_logger.info('Count information already collected in the output folder : %s --> Skipping this step!', rna_sc_count_path)
			df_counts = pd.read_csv(rna_sc_count_path)
			
			with open(alignment_information_sc_path, 'rb') as fp:
				alignment_information_sc = pickle.load(fp)

			df_alignments = pd.read_csv(rna_sc_count_all_alignments_path)

		try:
			os.remove(self.path_to_output_folder_alignments+'peptides_info.pkl')
		except FileNotFoundError:
			pass

		return df_counts, alignment_information_sc, df_alignments


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

		pos = alignment.split(':')[1].split('|')
		pos_set = []
		splice_pos = set()
		for i, chunck in enumerate(pos):
			ini = int(chunck.split('-')[0])
			fini = int(chunck.split('-')[1])
			aux = set(range(ini, fini+1))
			pos_set.extend(aux)
			if i != 0:
				splice_pos.add(ini)
			if len(pos) > 1 and i != len(pos) -1:
				splice_pos.add(fini)
			
		counts = self.get_depth_with_view(region_to_query, bam_file, name_sample, library, sens, strand, sequences, pos_set, splice_pos)
		

		to_return = [[peptide, alignment, name_sample, strand]]
		
		for index, sequence in enumerate(sequences):
			try:
				count = counts[sequence]
			except IndexError:
				count = -1
			to_return.append([count, sequence])
		
		return to_return


	def set_strand_read(self, strand):
		number = "{0:b}".format(int(strand))
		if len(number)>= 5:
			if '1' == number[-5]:
				return '-'
			else:
				return '+'
		else:
			return '+'

	def get_ranges(self, cigar, start, strand):
		rang = []
		splices_sites = set()
		rx = re.findall('(\d+)([MISDNX=])?', cigar)
		indels = []
		for index in rx:
			operation = index[1]
			length = int(index[0])
			
			if ('S' in operation):
				end = length
				start += end
				
			elif ('M' in index) or ('=' in index) or ('X' in index) :
				end = start + length
				rang.extend(range(start, end))
				start += length
				
			elif ('N' in index) or ('D' in index):
				splices_sites.add(end)
				start = start + length
				splices_sites.add(start)

			elif ('I' in index):
				indels.append(len(rang)+1)
		return rang, splices_sites, indels

	def get_indexes(self, overlap, rang_):
		return [rang_.index(i) for i in overlap]

	def remove_at(self, i, s):
		return s[:i] + s[i+1:]

	def define_read_in(self, read, pos_set, splice_pos):

		if len(read) > 0:
			split_read = read.split('\t')
			name = split_read[0]
			strand = split_read[1]
			start = int(split_read[3])
			seq = split_read[9]
			cigar = split_read[5]
			strand = self.set_strand_read(strand)
			rang_, splices_sites, indels = self.get_ranges(cigar, start, strand)
			splices_sites = splices_sites - splice_pos
			overlap = set(rang_).intersection(pos_set)

			try:
				cell = read.split('CR:Z:')[1].split('\t')[0]
			except:
				cell = read.split('CB:Z:')[1].split('-')[0]

			if len(indels) > 0:
				for indel in indels:
					seq = self.remove_at(indel, seq)

			if len(splices_sites.intersection(pos_set[1:-1])) == 0 and len(overlap) > 0:
				percentage_overlap = len(overlap)/len(pos_set)
				indexes = self.get_indexes(overlap, rang_)
				index_ini = min(indexes)
				index_fin = max(indexes) + 1
				seq_overlap = seq[index_ini :index_fin]
				if percentage_overlap == 1 :
					return name, cell, strand, seq_overlap, percentage_overlap
		
	def get_depth_with_view(self, region_to_query, bam_file, index_sample, library, sens, strand, sequences, pos_set, splice_pos):

		t_0 = time.time()
		contReads_to_return = {}

		library = library.lower()
		sens = sens.lower()
		
		def set_count(info_read, cs):

			name, cell, strand, seq_overlap, percentage_overlap = info_read
			
			if seq_overlap in cs and name not in set_names_reads:
				set_names_reads.add(name)
				cells_names_reads.add((name, cell))
				return percentage_overlap

		if library == 'unstranded' or sens == 'unstranded':

			try:
				count_1 = pysam.view("-F0X100", bam_file, region_to_query).split('\n')
			except pysam.utils.SamtoolsError: 
				return -1

			reads_overlaping_area = list(map(self.define_read_in, count_1, repeat(pos_set), repeat(splice_pos)))
			reads_overlaping_area = list(filter(None, reads_overlaping_area))
			reads_overlaping_area.sort(key=itemgetter(-1), reverse=True)

			
			for seq in sequences:
				contReads_to_return[seq] = {}
				contReads = 0
				set_names_reads = set()
				cells_names_reads = set()

				rcmcs = uf.reverseComplement(seq)
				
				sum_overlap_seq = list(map(set_count, reads_overlaping_area, repeat(seq)))
				
				for name, cell in cells_names_reads:
					name_cell = str(index_sample)+'_'+cell
					try:
						contReads_to_return[seq][name_cell] += 1
					except KeyError:
						contReads_to_return[seq][name_cell] = 1

				cells_names_reads = set()
				sum_overlap_seq = list(map(set_count, reads_overlaping_area, repeat(rcmcs)))
				
				for name, cell in cells_names_reads:
					name_cell = str(index_sample)+'_'+cell
					try:
						contReads_to_return[seq][name_cell] += 1
					except KeyError:
						contReads_to_return[seq][name_cell] = 1
				
		elif library == 'single-end':

			if ((strand == '+' and sens == 'forward') or (strand == '-' and sens == 'reverse')):
				try:
					count_1 = pysam.view("-F0X110", bam_file, region_to_query).split('\n')
				except pysam.utils.SamtoolsError: 
					return -1
			elif ((strand == '-' and sens == 'forward') or (strand == '+' and sens == 'reverse')):
				try:
					count_1 = pysam.view("-F0X100", "-f0X10", bam_file, region_to_query).split('\n')
				except pysam.utils.SamtoolsError: 
					return -1

			reads_overlaping_area = list(map(self.define_read_in, count_1, repeat(pos_set), repeat(splice_pos)))
			reads_overlaping_area = list(filter(None, reads_overlaping_area))
			reads_overlaping_area.sort(key=itemgetter(-1), reverse=True)

			for seq in sequences:
				if strand == '-':
					seq = uf.reverseComplement(seq)
				contReads_to_return[seq] = {}
				set_names_reads = set()
				cells_names_reads = set()
				sum_overlap_seq = list(map(set_count, reads_overlaping_area, repeat(seq)))
				
				for name, cell in cells_names_reads:
					name_cell = str(index_sample)+'_'+cell
					try:
						contReads_to_return[seq][name_cell] += 1
					except KeyError:
						contReads_to_return[seq][name_cell] = 1
				
		elif library == 'pair-end':
			if ((strand == '+' and sens == 'forward') or (strand == '-' and sens == 'reverse')):
				try:
					count_1 = pysam.view("-F0X100", "-f0X60", bam_file, region_to_query).split('\n')
					count_2 = pysam.view("-F0X100", "-f0X90", bam_file, region_to_query).split('\n')
				except pysam.utils.SamtoolsError: 
					return -1

			elif ((strand == '-' and sens == 'forward') or (strand == '+' and sens == 'reverse')):
				try:
					count_1 = pysam.view("-F0X100", "-f0X50", bam_file, region_to_query).split('\n')
					count_2 = pysam.view("-F0X100", "-f0XA0", bam_file, region_to_query).split('\n')
				except pysam.utils.SamtoolsError: 
					return -1

			reads_overlaping_area_count_1 = list(map(self.define_read_in, count_1, repeat(pos_set), repeat(splice_pos)))
			reads_overlaping_area_count_1 = list(filter(None, reads_overlaping_area_count_1))
			reads_overlaping_area_count_1.sort(key=itemgetter(-1), reverse=True)

			reads_overlaping_area_count_2 = list(map(self.define_read_in, count_2, repeat(pos_set), repeat(splice_pos)))
			reads_overlaping_area_count_2 = list(filter(None, reads_overlaping_area_count_2))
			reads_overlaping_area_count_2.sort(key=itemgetter(-1), reverse=True)

			for seq in sequences:
				contReads = 0
				set_names_reads = set()
				cells_names_reads = set()
				rcmcs_aux = uf.reverseComplement(seq)
				contReads_to_return[seq] = {}

				if strand == '-':
					seq = rcmcs_aux
					rcmcs = rcmcs_aux
				else:
					rcmcs = seq

				sum_overlap_seq = list(map(set_count, reads_overlaping_area_count_1, repeat(seq)))

				for name, cell in cells_names_reads:
					name_cell = str(index_sample)+'_'+cell
					try:
						contReads_to_return[seq][name_cell] += 1
					except KeyError:
						contReads_to_return[seq][name_cell] = 1

				cells_names_reads = set()
				sum_overlap_seq = list(map(set_count, reads_overlaping_area_count_2, repeat(rcmcs)))
				
				for name, cell in cells_names_reads:
					name_cell = str(index_sample)+'_'+cell
					try:
						contReads_to_return[seq][name_cell] += 1
					except KeyError:
						contReads_to_return[seq][name_cell] = 1
						
		# reading htslib https://github.com/DecodeGenetics/graphtyper/issues/57
		t_2 = time.time()
		total = t_2-t_0
		#print ('Total ', total)
		
		return contReads_to_return

