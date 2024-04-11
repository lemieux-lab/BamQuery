import itertools
import gc
import logging
import numpy as np
import os
import pandas as pd
import pickle
import pysam
import re
import sys
import time
import warnings
from concurrent.futures.process import _chain_from_iterable_of_lists
from pathos.multiprocessing import ProcessPool

import utils.useful_functions as uf
from genomics.get_counts_from_sample import get_counts_sample

if sys.version_info.minor >= 12:
	from itertools import batched

else:
	from itertools import islice

	def batched(iterable, n):
		if n < 1:
			raise ValueError('n must be at least one')
		it = iter(iterable)
		# starting from python v3.8, can use:
		# while batch := tuple(islice(it, n)):
		#	yield batch
		batch = tuple(itertools.islice(it, n))
		while batch:
			yield batch
			batch = tuple(islice(it, n))

warnings.filterwarnings("ignore")
os.environ["REF_PATH"] = "/dev/null"

__author__ = "Maria Virginia Ruiz Cuevas"
__email__ = "maria.virginia.ruiz.cuevas@umontreal.ca"


class GetCounts:

	def __init__(self, path_to_output_folder, name_exp, mode, light, peptides_by_type, super_logger, threads):
		self.light = light
		self.mode = mode
		self.path_to_output_folder = path_to_output_folder
		self.path_to_output_folder_alignments = path_to_output_folder+'alignments/'
		self.name_exp = name_exp
		self.peptides_by_type = peptides_by_type
		self.super_logger = super_logger
		self.threads = threads

		self.stdout_handler = logging.StreamHandler(sys.stdout)
		self.super_logger.addHandler(self.stdout_handler)

		if self.light:
			self.path_to_output_temps_folder = path_to_output_folder+'res_light/temps_files/'
			self.alignment_information_path = self.path_to_output_folder_alignments+'Alignments_information_light_rna.dic'

		else:
			self.path_to_output_temps_folder = path_to_output_folder+'res/temps_files/'
			self.alignment_information_path = self.path_to_output_folder_alignments+'Alignments_information_rna.dic'

		self.count_path = self.path_to_output_temps_folder+self.name_exp+'_rna_count.csv'
		self.count_all_alignments_path = self.path_to_output_temps_folder+self.name_exp+'_rna_count_All_alignments.csv'
		self.last_treated_bam_file = self.path_to_output_folder_alignments+'info_treated_bam_files.pkl'


		if self.mode == 'translation':
			self.path_to_output_temps_folder = path_to_output_folder+'res_translation/temps_files/'
			self.path_to_output_bed_folder = path_to_output_folder+'res_translation/BED_files/'
			self.count_path = self.path_to_output_temps_folder+self.name_exp+'_ribo_counts.csv'
			self.count_all_alignments_path = self.path_to_output_temps_folder+self.name_exp+'_ribo_count_All_alignments.csv'
			self.last_treated_bam_file = self.path_to_output_folder_alignments+'info_treated_ribo_bam_files.pkl'
			if self.light:
				self.alignment_information_path = self.path_to_output_folder_alignments+'Alignments_information_light_ribo.dic'
			else:
				self.alignment_information_path = self.path_to_output_folder_alignments+'Alignments_information_ribo.dic'


	def filter_alignment_information(self, perfect_alignments, path_alignment_information):

		alignment_information = {}
		for peptide_key, information_peptide in perfect_alignments.items():
			peptide = peptide_key.split('_')[0]
			if peptide in self.peptides_by_type :
				alignment_information[peptide_key] = information_peptide

		with open(path_alignment_information, 'wb') as handle:
			pickle.dump(alignment_information, handle, protocol=pickle.HIGHEST_PROTOCOL)

		return alignment_information


	def get_counts(self, perfect_alignments, bam_files_list, overlap):

		times = []
		df_counts = pd.DataFrame()

		exists_count = os.path.exists(self.count_path)
		exists_alignment_information = os.path.exists(self.alignment_information_path)

		last_treated_bam_file = os.path.exists(self.last_treated_bam_file)

		if last_treated_bam_file:

			with open(self.last_treated_bam_file , 'rb') as fp:
				last_treated_bam_file = pickle.load(fp)

		else:
			last_treated_bam_file = -1


		if not exists_count :
			t_0 = time.time()

			if not exists_alignment_information :
				alignment_information = self.filter_alignment_information(
					perfect_alignments,
					self.alignment_information_path)
			else:
				with open(self.alignment_information_path, 'rb') as fp:
			 		alignment_information = pickle.load(fp)

			info_bams = []
			bams = []

			list_bam_files_order = []
			for name_sample, info_bam in sorted(
					bam_files_list.items(),
					key=lambda e: e[1][-2],
					reverse=False
				):
				info_bams.append((name_sample,info_bam))
				list_bam_files_order.append(name_sample)

			def get_index(sample):
				index = list_bam_files_order.index(sample)
				return index + 1

			total_samples = len(list_bam_files_order)
			self.super_logger.info('Getting counts for %s samples', str(total_samples))

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

			infile = os.path.join(
				self.path_to_output_folder,
				'genome_alignments',
				'references_chrs.pkl'
			)
			with open(infile, "rb") as f:
				references_chrs = pickle.load(f)

			digits_bam_file_BQ_reference = []
			for chr in references_chrs:
				l = [x for x in chr if x.isdigit()]
				key = ''.join(l)
				if len(l) == 0:
					key = chr
				digits_bam_file_BQ_reference.append(key)

			if len(info_bams)-1 != last_treated_bam_file:
				pool = ProcessPool(nodes=self.threads)

			remaining_idxs = []
			remaining_bams = []

			for idx, bam_file in enumerate(info_bams):
				if idx > last_treated_bam_file:
					remaining_idxs.append(idx)
					remaining_bams.append(bam_file)
				else:
					self.super_logger.info(
						'Bam File already processed: %s %s.',
						str(idx), bam_file[0]
					)

			gc.collect()

			for idx, bam_file in zip(remaining_idxs, remaining_bams):
				if True:

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

					reference = self.get_corresponding_references(
						list(references),
						references_chrs,
						digits_bam_file_BQ_reference,
						bam_file_path
					)

					chunksize=100

					iresults = pool.imap(
						get_counts_sample,
						itertools.repeat(bam_file),
						batched(keys, chunksize),
						batched(values, chunksize),
						itertools.repeat(reference),
						itertools.repeat(overlap)
					)

					results = list(_chain_from_iterable_of_lists(iresults))

					message = "was successfully processed"

					for res in results:
						for index, count_align in enumerate(res):
							if index == 0:
								peptide = count_align[0]
								alignment = count_align[1]
								sample = count_align[2]
								strand = count_align[3]

							else:
								count = count_align[0]
								if count == -1:
									message = "couldn't be processed. Failed to open, permission denied"

								sequence = count_align[1]

								index_sample = get_index(sample)

								key = peptide+'_'+alignment+'_'+sequence
								if len(alignment_information[key][-1]) == 0:
									alignment_information[key][-1] = [0]*total_samples

								alignment_information[key][-1][index_sample-1] = count

					pool.clear()

					gc.collect()

					t1_bam_file = time.time()
					time_final = (t1_bam_file-t0_bam_file)/60.0
					times.append(time_final)

					self.super_logger.info(
						'Bam File : %s %s %s. Time : %f min',
						str(idx), bam_file[0], message, time_final
					)

					if (idx % 50 == 0) and (idx != 0):
						self.super_logger.info('Saving information for Bam Files processed')
						print(str(idx), 'Bam Files processed')

						with open(self.last_treated_bam_file, 'wb') as f:
							pickle.dump(idx, f)

						with open(self.alignment_information_path, 'wb') as handle:
							pickle.dump(alignment_information, handle, protocol=pickle.HIGHEST_PROTOCOL)

			if len(info_bams)-1 != last_treated_bam_file:
				pool.close()
				pool.join()
				pool.clear()

			gc.collect()

			self.super_logger.info('Average time to process a BamFile : %f min', np.mean(times))

			with open(self.last_treated_bam_file, 'wb') as f:
				pickle.dump(idx, f)

			with open(self.alignment_information_path, 'wb') as handle:
				pickle.dump(alignment_information, handle, protocol=pickle.HIGHEST_PROTOCOL)

			header = ['Peptide Type', 'Peptide', 'Position', 'MCS', 'Strand']
			header.extend(list_bam_files_order)
			to_write_list = []

			for peptide_key, info_peptide in alignment_information.items():

				split_key = peptide_key.split('_')
				peptide = split_key[0]
				position = split_key[1]
				MCS = split_key[2]
				strand = info_peptide[0]
				peptide_type = self.peptides_by_type[peptide]
				info_counts = info_peptide[-1]
				to_add =  [peptide_type, peptide, position, MCS, strand]
				to_add.extend(info_counts)
				to_write_list.append(to_add)

			df_alignments = pd.DataFrame(to_write_list, columns=header)
			df_counts = df_alignments.groupby(['Peptide Type', 'Peptide']).sum(numeric_only=True).reset_index()
			df_counts.sort_values(by=['Peptide Type'])
			df_counts.to_csv(self.count_path, header=True, index=False)

			df_alignments.to_csv(self.count_all_alignments_path, index=False, header=True)

			self.super_logger.info('Counts Information saved to : %s ', self.count_path)

			t_2 = time.time()
			total = t_2-t_0
			self.super_logger.info('Total time run function get_counts to end : %f min', (total/60.0))

		else:
			self.super_logger.info('Count information already collected in the output folder : %s --> Skipping this step!', self.count_path)
			df_counts = pd.read_csv(self.count_path)

			with open(self.alignment_information_path, 'rb') as fp:
				alignment_information = pickle.load(fp)

			df_alignments = pd.read_csv(self.count_all_alignments_path)

		self.super_logger.removeHandler(self.stdout_handler)
		self.stdout_handler.close()

		return df_counts, alignment_information, df_alignments


	def select_not_digits_no_special_chars(self, input_string):
		return "".join(re.findall(r'[^\d\W]', input_string))

	def is_contained(self, word, text):
		return bool(re.search(f'(?<=\w){word}(?=\w)', text))

	def has_digits(self, input_string):
		return bool(re.search(r'\d', input_string))

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
							if self.is_contained(self.select_not_digits_no_special_chars(references_chr), self.select_not_digits_no_special_chars(chr)):
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

