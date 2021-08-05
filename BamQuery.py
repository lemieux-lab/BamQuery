import time, sys, os, datetime, argparse, logging, _thread, shutil, pickle
from os import listdir
from os.path import isfile, join
import pandas as pd

from readers.read_input import ReadInputFile
from readers.info_transcripts_annotation import InfoTranscripts
from readers.intersection_alignments_annotations import IntersectAnnotations

from utils.primary_read_count import GetPrimaryReadCountBamFiles
from utils.reverse_translation import ReverseTranslation
from utils.paths_arrangements import *

from genomics.alignments import Alignments
from genomics.get_counts import GetCounts
from genomics.normalization import Normalization
from genomics.get_biotype import BiotypeAssignation

import plotting.plots as plots

__author__ = "Maria Virginia Ruiz Cuevas"
__email__ = "maria.virginia.ruiz.cuevas@umontreal.ca"


class BamQuery:

	def __init__(self, path_to_input_folder, path_to_output_folder, name_exp, mode, strandedness, th_out, light):
		self.path_to_input_folder = path_to_input_folder
		self.path_to_output_folder = path_to_output_folder
		self.name_exp = name_exp
		self.strandedness = strandedness
		self.mode = mode
		self.th_out = th_out
		self.light = light

		if self.mode == 'normal':
			self.run_bam_query_normal_mode()
		#else:
		#	self.run_bam_query_translation_mode()

		if not self.light:
			self.get_annotations()


	def run_bam_query_normal_mode(self):
		self.common_to_modes()
		
		if not self.light:
			name_path = self.path_to_output_folder+'/res/'+self.name_exp+'_count_norm_info.xlsx'
		else:
			name_path = self.path_to_output_folder+'/res_light/'+self.name_exp+'_count_norm_info.xlsx'
		
		exists = os.path.exists(name_path) 

		if not exists:
			
			get_counts = GetCounts(self.path_to_output_folder, self.name_exp, self.mode, self.light, self.input_file_treatment.peptides_by_type)
			
			res = get_counts.get_counts(self.perfect_alignments, self.bam_files_info.bam_files_list)
			df_counts_rna = res[0] 
			self.perfect_alignments = res[1]
			df_counts_filtered = res[2] 
			order = res[3] 
			order_f = res[4] 
			df_all_alignments_rna = res[5] 

			writer = pd.ExcelWriter(name_path, engine='xlsxwriter')
			df_all_alignments_rna.to_excel(writer, sheet_name='Alignments Read count RNA-seq')
			df_counts_rna.to_excel(writer, sheet_name='Read count RNA-seq by peptide')
			
			if not self.light: 
				plots.get_heat_map(df_counts_rna, self.path_to_output_folder, self.name_exp, '_rna_counts', False, order, self.th_out)

			logging.info('========== Get Count RNA : Done! ============ ')

			normalization = Normalization(self.path_to_output_folder, self.name_exp, self.bam_files_info.bam_files_list, self.input_file_treatment.all_mode_peptide, self.mode, self.light)
			def_norm_rna = normalization.get_normalization(df_counts_rna, '_rna_norm.csv')
			def_norm_rna.to_excel(writer, sheet_name='log10(RPHM) RNA-seq by peptide')
			
			if not self.light: 
				plots.get_heat_map(def_norm_rna, self.path_to_output_folder, self.name_exp, '_rna_norm', True, order, self.th_out)

			writer.save()
			logging.info('========== Get Norm RNA : Done! ============ ')
		else:
			print ('Information count and normalisation already collected !')
		
	def run_bam_query_translation_mode(self):
		self.common_to_modes()

		exists = os.path.exists(self.path_to_output_folder+'/res/'+self.name_exp+'_count_norm_info.xlsx')

		writer = pd.ExcelWriter(self.path_to_output_folder+'/res/'+self.name_exp+'_count_norm_info.xlsx', engine='xlsxwriter')

		get_counts = GetCounts(self.path_to_output_folder, self.name_exp, self.mode, self.input_file_treatment.peptides_by_type)
		perfect_alignments_to_return, df_counts_ribo, order, df_all_alignments_ribo = get_counts.ribo_counts(self.perfect_alignments, self.bam_files_info.bam_ribo_files_list)
		df_all_alignments_ribo.to_excel(writer, sheet_name='Alignments Read count Ribo-seq')
		df_counts_ribo.to_excel(writer, sheet_name='Read count Ribo-seq by peptide')

		if len(df_counts_ribo) < 400:
			plots.get_heat_map(df_counts_ribo, self.path_to_output_folder, self.name_exp, '_ribo_counts', False, order)
		
		logging.info('========== Get Count Ribo : Done! ============ ')
		print ('Get Count Ribo : Done!')

		normalization = Normalization(self.path_to_output_folder, self.name_exp, self.bam_files_info.bam_ribo_files_list, self.input_file_treatment.all_mode_peptide, self.mode)
		def_norm_ribo = normalization.get_normalization(df_counts_ribo, '_ribo_norm.csv')
		
		def_norm_ribo.to_excel(writer, sheet_name='log10(RPHM) Ribo-seq by peptide')
		if len(def_norm_ribo) < 400:
			plots.get_heat_map(def_norm_ribo, self.path_to_output_folder, self.name_exp, '_ribo_norm', True, order, self.th_out)

		logging.info('========== Get Norm Ribo : Done! ============ ')
		print ('Get Norm Ribo : Done!')

		df_counts_rna, self.perfect_alignments, df_counts_filtered, order, order_f, df_all_alignments_rna = get_counts.get_counts(perfect_alignments_to_return, self.bam_files_info.bam_files_list)
		df_all_alignments_rna.to_excel(writer, sheet_name='Alignments Read count RNA-seq')
		df_counts_rna.to_excel(writer, sheet_name='Read count RNA-seq by peptide')

		if len(df_counts_rna) < 400:
			plots.get_heat_map(df_counts_rna, self.path_to_output_folder, self.name_exp, '_rna_counts', False, order)
		
		if len(df_counts_filtered) < 400:
			plots.get_heat_map(df_counts_filtered, self.path_to_output_folder, self.name_exp, '_rna_ribo_counts', False, order_f)

		logging.info('========== Get Count RNA : Done! ============ ')
		print ('Get Count RNA : Done!')

		normalization = Normalization(self.path_to_output_folder, self.name_exp, self.bam_files_info.bam_files_list, self.input_file_treatment.all_mode_peptide, self.mode)
		def_norm_rna = normalization.get_normalization(df_counts_rna, '_rna_norm.csv')
		def_norm_rna.to_excel(writer, sheet_name='log10(RPHM) RNA-seq by peptide')
		plots.get_heat_map(def_norm_rna, self.path_to_output_folder, self.name_exp, '_rna_norm', True, order, self.th_out)

		logging.info('========== Get Norm RNA : Done! ============ ')
		print ('Get Norm RNA : Done!')

		def_norm_rna_ribo = normalization.get_normalization(df_counts_filtered, '_rna_ribo_norm.csv')
		def_norm_rna_ribo.to_excel(writer, sheet_name='log10(RPHM) RNA and Ribo counts')
		
		if len(def_norm_rna_ribo) < 400:
			plots.get_heat_map(def_norm_rna_ribo, self.path_to_output_folder, self.name_exp, '_rna_ribo_norm', True, order_f, self.th_out)

		logging.info('========== Get Norm Ribo-RNA : Done! ============ ')
		print ('Get Norm Ribo-RNA : Done!')

		writer.save()
		

	def common_to_modes(self):
		self.bam_files_info = GetPrimaryReadCountBamFiles()
		self.bam_files_info.get_all_counts(self.path_to_input_folder, self.path_to_output_folder, self.mode, self.strandedness, self.light)

		logging.info('========== Get Primary Counts : Done! ============ ')
		print ('Get Primary Counts : Done!')

		self.input_file_treatment = ReadInputFile(self.path_to_input_folder)
		self.input_file_treatment.treatment_file()


		with open(self.path_to_output_folder+'genome_alignments/peptides_by_type.dic', 'wb') as handle:
			pickle.dump(self.input_file_treatment.peptides_by_type, handle, protocol=pickle.HIGHEST_PROTOCOL)

		with open(self.path_to_output_folder+'genome_alignments/all_mode_peptide.dic', 'wb') as handle:
			pickle.dump(self.input_file_treatment.all_mode_peptide, handle, protocol=pickle.HIGHEST_PROTOCOL)

		
		logging.info('========== Treatment File : Done! ============ ')
		print ('Treatment File : Done!')
		set_peptides = set(list(self.input_file_treatment.all_mode_peptide.keys()))

		if not self.light:
			name_path = self.path_to_output_folder+'alignments/Alignments_information.dic'
		else:
			name_path = self.path_to_output_folder+'alignments/Alignments_information_light.dic'
		
		exists = os.path.exists(name_path) 

		if not exists:
			if len(self.input_file_treatment.peptide_mode) > 0 or len(self.input_file_treatment.CS_mode) > 0 :

				if len(self.input_file_treatment.CS_mode) > 0 :			
					with open(self.path_to_output_folder+'genome_alignments/CS_mode.dic', 'wb') as handle:
						pickle.dump(self.input_file_treatment.CS_mode, handle, protocol=pickle.HIGHEST_PROTOCOL)

				if len(self.input_file_treatment.peptide_mode) > 0 :			
					with open(self.path_to_output_folder+'genome_alignments/peptides_mode.dic', 'wb') as handle:
						pickle.dump(self.input_file_treatment.peptide_mode, handle, protocol=pickle.HIGHEST_PROTOCOL)

				self.reverse_translation = ReverseTranslation()
				self.reverse_translation.reverse_translation(self.input_file_treatment.peptide_mode, self.input_file_treatment.CS_mode, self.path_to_output_folder, self.name_exp)
				
				logging.info('========== Reverse Translation : Done! ============ ')
				print ('Reverse Translation : Done!')
			
				self.alignments = Alignments(self.path_to_output_folder, self.name_exp, self.light)
				self.perfect_alignments, peptides_with_alignments = self.alignments.alignment_cs_to_genome(set_peptides)

				logging.info('========== Alignment : Done! ============ ')
				print ('Alignment : Done!')

			else:
				self.perfect_alignments = {}
				peptides_with_alignments = set()

				if not self.light:
					name_path = self.path_to_output_folder+'alignments/Alignments_information.dic'
				else:
					name_path = self.path_to_output_folder+'alignments/Alignments_information_light.dic'

				with open(name_path, 'wb') as handle:
					pickle.dump(self.perfect_alignments, handle, protocol=pickle.HIGHEST_PROTOCOL)
		else:
			with open(name_path, 'rb') as handle:
				self.perfect_alignments = pickle.load(handle)


		# positions_mcs_peptides_variants_alignment[key] = [strand, local_translation_peptide, differences_pep, info_snps, differences_ntds, [],[]]
		exists = os.path.exists(self.path_to_output_folder+'/res/'+self.name_exp+'_count_norm_info.xlsx') 

		if not exists:
			if len(self.input_file_treatment.manual_mode) > 0 :
				with open(self.path_to_output_folder+'genome_alignments/manual_mode.dic', 'wb') as handle:
					pickle.dump(self.input_file_treatment.manual_mode, handle, protocol=pickle.HIGHEST_PROTOCOL)

				for peptide, info_peptide in self.input_file_treatment.manual_mode.items() :
					coding_sequence = info_peptide[0]
					position = info_peptide[1]
					strand = info_peptide[2]
					key = peptide+'_'+position+'_'+coding_sequence
					peptides_with_alignments.add(peptide)
					self.perfect_alignments[key] = [strand, peptide, ['NA'], ['NA'], ['NA'], [], []]

		exists = os.path.exists(self.path_to_output_folder+'alignments/missed_peptides.info')

		if not exists:
			missed_peptides = list(set_peptides - peptides_with_alignments)

			with open(self.path_to_output_folder+'alignments/missed_peptides.info', 'w') as f:
				for item in missed_peptides:
					f.write("%s\t\n" % item)

			logging.info('Total missed_peptides : %s. Find the list in : %s.', str(len(missed_peptides)), self.path_to_output_folder+'alignments/missed_peptides.info')

		logging.info('========== Common_to_modes : Done! ============ ')
		print ('common_to_modes : Done!')


	def get_annotations(self):

		logging.info('========== Running get_annotations ============ ')
		info_peptide_alignments = self.get_info_peptide_alignments()

		get_info_transcripts = InfoTranscripts()
		get_info_transcripts.set_values()

		if len(self.perfect_alignments) == 0:
			if not self.light:
				name_path = self.path_to_output_folder+'alignments/Alignments_information.dic'
			else:
				name_path = self.path_to_output_folder+'alignments/Alignments_information_light.dic'

			with open(name_path, 'rb') as fp:
				self.perfect_alignments = pickle.load(fp)

		intersect_to_annotations = IntersectAnnotations(self.perfect_alignments, self.path_to_output_folder, self.name_exp)
		intersect_to_annotations.generate_BED_files()
		intersect_to_annotations.perform_intersection_with_annotation()

		logging.info('========== Intersect to Annotations : Done! ============ ')

		list_bam_files_order_rna = []
		order_sample_bam_files_rna = {}

		for name_sample, info_bam in sorted(self.bam_files_info.bam_files_list.items(), key=lambda e: e[1][-1], reverse=False):
			list_bam_files_order_rna.append(name_sample)
			group = info_bam[-1]
			try:
				order_sample_bam_files_rna[group].append(name_sample)
			except KeyError:
				order_sample_bam_files_rna[group] = [name_sample]

		list_bam_files_order_ribo = []
		order_sample_bam_files_ribo = {}

		for name_sample, info_bam in sorted(self.bam_files_info.bam_ribo_files_list.items(), key=lambda e: e[1][-1], reverse=False):
			list_bam_files_order_ribo.append(name_sample)
			try:
				order_sample_bam_files_ribo[group].append(name_sample)
			except KeyError:
				order_sample_bam_files_ribo[group] = [name_sample]
		
		get_biotype = BiotypeAssignation(self.path_to_output_folder, self.name_exp, self.mode, list_bam_files_order_rna, list_bam_files_order_ribo, order_sample_bam_files_rna, order_sample_bam_files_ribo)
		
		get_biotype.get_biotypes(info_peptide_alignments, self.input_file_treatment.peptides_by_type)
		get_biotype.get_global_annotation()
		
		logging.info('========== Annotations : Done! ============ ')


	def get_info_peptide_alignments(self):

		info_peptide_alignments = {}

		for peptide_alignment in self.perfect_alignments:
			peptide = peptide_alignment.split('_')[0]
			alignment = peptide_alignment.split('_')[1]
			MCS = peptide_alignment.split('_')[2]
			count_rna = self.perfect_alignments[peptide_alignment][-2]
			count_ribo = self.perfect_alignments[peptide_alignment][-1]
			strand = self.perfect_alignments[peptide_alignment][0]
			
			try:
				info_peptide = info_peptide_alignments[peptide]
				info_peptide[0][peptide_alignment] = [strand, count_rna, count_ribo]
				info_peptide[1] += sum(count_rna) 
				info_peptide[2] += sum(count_ribo)
			except KeyError:
				dic = {}
				dic[peptide_alignment] = [strand, count_rna, count_ribo]
				info_peptide_alignments[peptide] = [dic, sum(count_rna), sum(count_ribo)]

		return info_peptide_alignments


def main(argv):

	parser = argparse.ArgumentParser(description='======== BamQuery ========')
	
	parser.add_argument('path_to_input_folder', type=str,
						help='Path to the input folder where to find BAM_directories.tsv and peptides.tsv')
	parser.add_argument('name_exp', type=str,
						help='BamQuery search Id')
	parser.add_argument('--mode', type=str, default = 'normal',
						help='BamQuery search mode : normal / translation')
	parser.add_argument('--strandedness', action='store_true',
						help='Take into account strandedness of the samples')
	parser.add_argument('--th_out', type=float, default = 8.55,
						help='Threshold to assess expression comparation with other tissues')
	parser.add_argument('--light', action='store_true',
						help='Display only the count and norm count for peptides and regions')

	args = parser.parse_args()
	
	path_to_input_folder = args.path_to_input_folder
	name_exp = args.name_exp.lower()
	mode = args.mode.lower()
	strandedness = args.strandedness
	th_out = args.th_out
	light = args.light

	if path_to_input_folder[-1] != '/':
		path_to_input_folder += '/'

	path_to_output_folder = directories_creation(path_to_input_folder, name_exp, mode, strandedness, light)

	t0 = time.time()

	BamQuery(path_to_input_folder, path_to_output_folder, name_exp, mode, strandedness, th_out, light)
	
	#try:
	#	shutil.rmtree(path_to_output_folder+'genome_alignments')
	#except FileNotFoundError:
	#	pass

	# if not light:
	# 	try:
	# 		shutil.rmtree(path_to_output_folder+'res/BED_files')
	# 		shutil.rmtree(path_to_output_folder+'res/AUX_files')
	# 		shutil.rmtree(path_to_output_folder+'res/temps_files')
	# 	except FileNotFoundError:
	# 		pass
	# else:
	# 	try:
	# 		shutil.rmtree(path_to_output_folder+'res_light/AUX_files')
	# 		shutil.rmtree(path_to_output_folder+'res_light/temps_files')
	# 		shutil.rmtree(path_to_output_folder+'res_light/BED_files')
	# 	except FileNotFoundError:
	# 		pass

	t2 = time.time()
	total = t2-t0
	
	logging.info('Total time run function BamQuery to end : %f min', (total/60.0))


if __name__ == "__main__":
	main(sys.argv[1:])


