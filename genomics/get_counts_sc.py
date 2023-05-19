import warnings, os
warnings.filterwarnings("ignore")
import os, time, pickle, os, pysam
import pandas as pd
from pathos.multiprocessing import ProcessPool
import utils.useful_functions as uf
from itertools import repeat
from operator import itemgetter
import re
from functools import reduce
import sys

__author__ = "Maria Virginia Ruiz Cuevas"
__email__ = "maria.virginia.ruiz.cuevas@umontreal.ca"


class GetCountsSC:

	def __init__(self, path_to_output_folder, name_exp, mode, light, peptides_by_type, super_logger, threads):
		self.path_to_output_folder = path_to_output_folder
		self.path_to_output_folder_res = path_to_output_folder+'res/'
		self.name_exp = name_exp
		self.mode = mode
		self.path_to_output_folder_alignments = path_to_output_folder+'alignments/'
		self.peptides_by_type = peptides_by_type
		self.super_logger = super_logger
		self.threads = threads

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

		df_counts = pd.DataFrame()
		
		rna_sc_count_path = self.path_to_output_folder_res+self.name_exp+'_rna_sc_count.csv'
		rna_sc_count_all_alignments_path = self.path_to_output_folder_res+self.name_exp+'_rna_sc_count_All_alignments.csv'
		
		exists_rna_sc = os.path.exists(rna_sc_count_path)
		
		alignment_information_sc_path = self.path_to_output_folder_alignments+'Alignments_information_sc.dic'
		exists_alignment_information_sc = os.path.exists(alignment_information_sc_path)
		
		last_treated_bam_file = os.path.exists(self.path_to_output_folder_alignments+'info_treated_bam_files.pkl')
		
		if last_treated_bam_file:

			with open(self.path_to_output_folder_alignments+'info_treated_bam_files.pkl', 'rb') as fp:
				last_treated_bam_file = pickle.load(fp)

			try:
				with open(self.path_to_output_folder_alignments+'peptides_info.pkl', 'rb') as fp:
					peptides_info = pickle.load(fp)

			except FileNotFoundError:
				peptides_info = {}

			try:
				with open(self.path_to_output_folder_alignments+'cell_names.pkl', 'rb') as fp:
					cell_names = pickle.load(fp)

			except FileNotFoundError:
				cell_names = set()

		else:
			last_treated_bam_file = -1
			peptides_info = {}
			cell_names = set()

		if not exists_rna_sc:
			t_0 = time.time()
			
			if not exists_alignment_information_sc :
				alignment_information_sc = self.filter_alignment_information(perfect_alignments, alignment_information_sc_path)
			else:
				with open(alignment_information_sc_path, 'rb') as fp:
			 		alignment_information_sc = pickle.load(fp)

			info_bams = []
			info_bams_path = []
			bams = []
			
			list_bam_files_order = []
			
			index_sample = 0
			for name_sample, info_bam in sorted(bam_files_list.items(), key=lambda e: e[1][-2], reverse=False):
				info_bams.append((index_sample,info_bam))
				info_bams_path.append((index_sample, name_sample, info_bam[0]))
				list_bam_files_order.append(name_sample)
				index_sample += 1
				
			df_info_bams = pd.DataFrame(info_bams_path, columns=['Index','Name Sample', 'Bam path'])
			df_info_bams.to_csv(self.path_to_output_folder_res+self.name_exp+'_index_bam_files.csv', index=None)
			
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

			pool = ProcessPool(nodes=self.threads)

			with open(self.path_to_output_folder+'/genome_alignments/references_chrs.pkl', "rb") as f:
				references_chrs = pickle.load(f)

			digits_bam_file_BQ_reference = []
			for chr in references_chrs:
				l = [x for x in chr if x.isdigit()]
				key = ''.join(l)
				if len(l) == 0:
					key = chr
				digits_bam_file_BQ_reference.append(key)

			for idx, bam_file in enumerate(info_bams):
				
				if idx > last_treated_bam_file:

					t0_bam_file = time.time()

					bam_file_path = bam_file[1][0]
					try:
						bam_file_ref = pysam.AlignmentFile(bam_file_path, "rb")
						references = bam_file_ref.header.references
						bam_file_ref.close()
					except ValueError as e:
						bam_file_ref = pysam.AlignmentFile(bam_file_path, "rb", check_sq=False)
						references = set()
						for alignment in bam_file_ref.fetch():
							references.add(bam_file_ref.getrname(alignment.reference_id))
						bam_file_ref.close()

					reference = self.get_corresponding_references(list(references), references_chrs, digits_bam_file_BQ_reference, bam_file_path)
					references = [reference] * len(keys)
					bams = [bam_file] * len(keys)
					results = pool.map(self.get_counts_sample, bams, keys, values, references)
					not_permission = False

					for res in results:
						for index, count_align in enumerate(res):

							if index == 0: 
								peptide = count_align[0]
								alignment = count_align[1]
								index_sample = count_align[2]
								strand = count_align[3]
								key = alignment+'_'+strand
							else:
								count_info = count_align[0]
								sequence = count_align[1]
								new_key = peptide+'_'+alignment+'_'+strand+'_'+sequence
								
								if count_info == -1:
									not_permission = True
								
								if len(count_info) > 0: # cell names
									cells = set(count_info.keys())
									cell_names = cell_names.union(cells)
								
								try:
									info_sequence = peptides_info[new_key]

									if len(count_info) > 0:
										info_sequence.update(count_info)
										
								except KeyError:
									peptides_info[new_key] = count_info
									
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

						with open(self.path_to_output_folder_alignments+'cell_names.pkl', 'wb') as f:  
							pickle.dump(cell_names, f)

						with open(alignment_information_sc_path, 'wb') as handle:
							pickle.dump(alignment_information_sc, handle, protocol=pickle.HIGHEST_PROTOCOL)

						with open(self.path_to_output_folder_alignments+'peptides_info.pkl', 'wb') as handle:
							pickle.dump(peptides_info, handle, protocol=pickle.HIGHEST_PROTOCOL)
				else:
					self.super_logger.info('Bam File already processed: %s %s.', str(idx), list_bam_files_order[idx])
			
			pool.close()
			pool.join()
			pool.clear()

			with open(self.path_to_output_folder_alignments+'info_treated_bam_files.pkl', 'wb') as f:  
				pickle.dump(idx, f)
			
			with open(self.path_to_output_folder_alignments+'cell_names.pkl', 'wb') as f:  
				pickle.dump(cell_names, f)

			with open(alignment_information_sc_path, 'wb') as handle:
				pickle.dump(alignment_information_sc, handle, protocol=pickle.HIGHEST_PROTOCOL)

			with open(self.path_to_output_folder_alignments+'peptides_info.pkl', 'wb') as handle:
				pickle.dump(peptides_info, handle, protocol=pickle.HIGHEST_PROTOCOL)

			cell_names = list(cell_names)
			if len(cell_names) == 0:
				header = ['Peptide Type', 'Peptide', 'Position', 'MCS', 'Strand', 'No intersected cells']
			else:
				cell_names = sorted(cell_names, key=lambda x: (int(x.split('_')[0]), x))
				data = []
				new_cell_names = []

				for cont, cell_name in enumerate(cell_names, start=1):
					bam_file = int(cell_name.split('_')[0])
					cell_name = cell_name.split('_')[1]
					new_cell_names.append(cont)
					data.append([cont, bam_file, cell_name])

				convertion = pd.DataFrame(data, columns=['id','Single_cell_bam_file', 'Cell name'])
				convertion.to_csv(self.path_to_output_folder_res+self.name_exp+'_cell_names_identification.csv',  header=True, index=False)
				header = ['Peptide Type', 'Peptide', 'Position', 'MCS', 'Strand']+new_cell_names
			
			data = []
			
			for peptide_info, cell_counts in peptides_info.items():
				peptide_info_split = peptide_info.split('_')
				peptide = peptide_info_split[0]
				alignment = peptide_info_split[1]
				strand = peptide_info_split[2]
				sequence = peptide_info_split[3]
				peptide_type = self.peptides_by_type[peptide]
				if len(cell_counts)>0 :  # remove the = if only want to show the peptides that have some counting
					aux = [peptide_type, peptide, alignment, sequence, strand]
					zeros = [0]*len(cell_names)
					for cell, count in cell_counts.items():
						index = cell_names.index(cell)
						zeros[index] = count
					aux = aux+zeros
					data.append(aux)
				if len(cell_names) == 0:
					aux = [peptide_type, peptide, alignment, sequence, strand, 0]
					data.append(aux)

			try:
				df_alignments = pd.DataFrame(data, columns = header)
			except MemoryError:
				print ('The available memory is insufficient to complete the query. Please allocate additional memory to ensure successful execution. \nIf the issue continues, we recommend dividing the number of peptides into smaller batches for the query.')
				sys.exit(2)

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

		return df_counts, alignment_information_sc, df_alignments


	def get_counts_sample(self, bam, peptide_alignment, sequences, references):
		
		name_sample = bam[0]
		bam_file = bam[1][0]
		library = bam[1][1]
		sens = bam[1][2]

		count = 0
		peptide = peptide_alignment.split('_')[0]
		alignment = peptide_alignment.split('_')[1]
		strand = peptide_alignment.split('_')[2]
		to_return = [[peptide, alignment, name_sample, strand]]

		chr = alignment.split(':')[0]
		try:
			chr = references[chr]
		except:
			for sequence in sequences:
				to_return.append([0, sequence])
			return to_return
		
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

	def get_ranges(self, cigar, start, len_seq):
		rang = [0]*len_seq
		splices_sites = set()
		rx = re.findall('(\d+)([MISDNX=])?', cigar)
		indels = []
		lastIndex = 0

		for index in rx:
			operation = index[1]
			length = int(index[0])
			
			if ('S' in operation):
				end = length
				lastIndex += end
				
			elif ('M' in index) or ('=' in index) or ('X' in index) :
				end = length

				for i in range(0,end):
					rang[lastIndex+i] = start
					start = start + 1

				lastIndex = lastIndex + end
				
			elif ('N' in index) or ('D' in index):
				splices_sites.add(start)
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
			rang_, splices_sites, indels = self.get_ranges(cigar, start, len(seq))
			splices_sites = splices_sites - splice_pos
			overlap = set(rang_).intersection(pos_set)
			percentage_overlap = len(overlap)/len(pos_set)
			
			try:
				cell = read.split('CB:Z:')[1].split('-')[0]
			except:
				cell = read.split('CR:Z:')[1].split('\t')[0]

			if len(indels) > 0:
				for indel in indels:
					seq = self.remove_at(indel, seq)

			if len(splices_sites.intersection(sorted(pos_set)[1:-1])) == 0 and percentage_overlap >= 0.6:
				index_ini = rang_.index(min(overlap))
				index_fin = rang_.index(max(overlap)) + 1
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
		
		return contReads_to_return
	
	def get_corresponding_references(self, references, references_chrs, digits_bam_file_BQ_reference, bam_file_path):
		reference_correspondence = {}
		digits_bam_file_search = {}
		no_digits = []
		chr_string = 'chr' in references[0]
		not_associated_chr_ref_bam_file = []
		
		for chr in references:
			l = [x for x in chr if x.isdigit()]
			key = ''.join(l)
			if len(l) == 0:
				key = chr
				no_digits.append(chr)
			try:
				digits_bam_file_search[key].append(chr)
			except KeyError:
				digits_bam_file_search[key] = [chr]
				
		for index, digits_reference in enumerate(digits_bam_file_BQ_reference):
			references_chr = references_chrs[index]
			try:
				chr_bam_file = digits_bam_file_search[digits_reference]
				if len(chr_bam_file) > 1:
					intersection = set([references_chr]).intersection(set(chr_bam_file))
					if intersection:
						reference_correspondence[references_chr] = list(intersection)[0]
					else:
						for chr in chr_bam_file:
							if self.is_contained_ignore_case(self.select_not_digits_no_special_chars(references_chr), self.select_not_digits_no_special_chars(chr)):
								reference_correspondence[references_chr] = chr
								break
						not_associated_chr_ref_bam_file.append(references_chr)
				else:
					reference_correspondence[references_chr] = chr_bam_file[0]
				
			except:
				in_ = False
				if not self.has_digits(digits_reference):
					if not chr_string:
						for index, i in enumerate(no_digits):
							if i == 'MT' and references_chr == 'chrM':
								reference_correspondence[references_chr] = i
								in_ = True
								no_digits.pop(index)
							else:
								chr_name = 'chr'+i
								if chr_name == references_chr:
									reference_correspondence[references_chr] = i
									in_ = True
									no_digits.pop(index)
							if in_:
								break
								
					else:
						for index, i in enumerate(no_digits):
							if ('MT' in i or 'M' in i) and references_chr == 'chrM':
								reference_correspondence[references_chr] = i
								in_ = True
								no_digits.pop(index)
							else: 
								if i.lower() == references_chr.lower():
									reference_correspondence[references_chr] = i
									in_ = True
									no_digits.pop(index)
							if in_:
								break

				if not in_:
					for digits_bam_file, list_chrs in digits_bam_file_search.items():
						if digits_reference in digits_bam_file:
							reference_correspondence[references_chr] = list_chrs[0]
							in_ = True
							break
				if not in_:
					not_associated_chr_ref_bam_file.append(references_chr)

		if len(not_associated_chr_ref_bam_file):
			chrs_not_asso = ','.join(not_associated_chr_ref_bam_file)
			self.super_logger.info('It was not possible to associate these chromosomes in the genome annotations reference: %s to any chr in the bam file investigated: %s. \nRegions in these chromosomes will be associated with a count of zero.', chrs_not_asso, bam_file_path)

		return reference_correspondence

