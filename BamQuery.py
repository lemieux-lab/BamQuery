import time, sys, os, datetime, argparse, logging, _thread
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

	def __init__(self, path_to_input_folder, path_to_output_folder, name_exp, mode, strandedness, th_out):
		self.path_to_input_folder = path_to_input_folder
		self.path_to_output_folder = path_to_output_folder
		self.name_exp = name_exp
		self.strandedness = strandedness
		self.mode = mode
		self.th_out = th_out

		if self.mode == 'normal':
			self.run_bam_query_normal_mode()
		else:
			self.run_bam_query_translation_mode()

		self.get_annotations()


	def run_bam_query_normal_mode(self):
		self.common_to_modes()
		
		writer = pd.ExcelWriter(self.path_to_output_folder+'/res/'+self.name_exp+'_count_norm_info.xlsx', engine='xlsxwriter')

		get_counts = GetCounts(self.path_to_output_folder, self.name_exp, self.mode)
		df_counts_rna, self.perfect_alignments, df_counts_filtered, df_all_alignments_rna = get_counts.get_counts(self.perfect_alignments, self.bam_files_info.bam_files_list)
		df_all_alignments_rna.to_excel(writer, sheet_name='Alignments Read count RNA-seq')
		df_counts_rna.to_excel(writer, sheet_name='Read count RNA-seq by peptide')

		plots.get_heat_map(df_counts_rna, self.path_to_output_folder, self.name_exp, '_rna_counts', False)

		normalization = Normalization(self.path_to_output_folder, self.name_exp, self.bam_files_info.bam_files_list, self.input_file_treatment.all_mode_peptide, self.mode)
		def_norm_rna = normalization.get_normalization(df_counts_rna, '_rna_norm.csv')
		def_norm_rna.to_excel(writer, sheet_name='log10(RPHM) RNA-seq by peptide')
		plots.get_heat_map(def_norm_rna, self.path_to_output_folder, self.name_exp, '_rna_norm', True, self.th_out)

		writer.save()
	
	def run_bam_query_translation_mode(self):
		self.common_to_modes()
		writer = pd.ExcelWriter(self.path_to_output_folder+'/res/'+self.name_exp+'_count_norm_info.xlsx', engine='xlsxwriter')
		

		get_counts = GetCounts(self.path_to_output_folder, self.name_exp, self.mode)
		perfect_alignments_to_return, df_counts_ribo, df_all_alignments_ribo = get_counts.ribo_counts(self.perfect_alignments, self.bam_files_info.bam_ribo_files_list)
		df_all_alignments_ribo.to_excel(writer, sheet_name='Alignments Read count Ribo-seq')
		df_counts_ribo.to_excel(writer, sheet_name='Read count Ribo-seq by peptide')

		plots.get_heat_map(df_counts_ribo, self.path_to_output_folder, self.name_exp, '_ribo_counts', False)
		
		logging.info('========== Get Count Ribo : Done! ============ ')
		print ('Get Count Ribo : Done!')

		normalization = Normalization(self.path_to_output_folder, self.name_exp, self.bam_files_info.bam_ribo_files_list, self.input_file_treatment.all_mode_peptide, self.mode)
		def_norm_ribo = normalization.get_normalization(df_counts_ribo, '_ribo_norm.csv')
		
		def_norm_ribo.to_excel(writer, sheet_name='log10(RPHM) Ribo-seq by peptide')
		
		plots.get_heat_map(def_norm_ribo, self.path_to_output_folder, self.name_exp, '_ribo_norm', True, self.th_out)

		logging.info('========== Get Norm Ribo : Done! ============ ')
		print ('Get Norm Ribo : Done!')

		df_counts_rna, self.perfect_alignments, df_counts_filtered, df_all_alignments_rna = get_counts.get_counts(perfect_alignments_to_return, self.bam_files_info.bam_files_list)
		df_all_alignments_rna.to_excel(writer, sheet_name='Alignments Read count RNA-seq')
		df_counts_rna.to_excel(writer, sheet_name='Read count RNA-seq by peptide')

		plots.get_heat_map(df_counts_rna, self.path_to_output_folder, self.name_exp, '_rna_counts', False)
		
		plots.get_heat_map(df_counts_filtered, self.path_to_output_folder, self.name_exp, '_rna_ribo_counts', False)

		logging.info('========== Get Count RNA : Done! ============ ')
		print ('Get Count RNA : Done!')

		normalization = Normalization(self.path_to_output_folder, self.name_exp, self.bam_files_info.bam_files_list, self.input_file_treatment.all_mode_peptide, self.mode)
		def_norm_rna = normalization.get_normalization(df_counts_rna, '_rna_norm.csv')
		def_norm_rna.to_excel(writer, sheet_name='log10(RPHM) RNA-seq by peptide')
		plots.get_heat_map(def_norm_rna, self.path_to_output_folder, self.name_exp, '_rna_norm', True, self.th_out)

		logging.info('========== Get Norm RNA : Done! ============ ')
		print ('Get Norm RNA : Done!')

		def_norm_rna_ribo = normalization.get_normalization(df_counts_filtered, '_rna_ribo_norm.csv')
		def_norm_rna_ribo.to_excel(writer, sheet_name='log10(RPHM) RNA and Ribo counts')
		
		plots.get_heat_map(def_norm_rna_ribo, self.path_to_output_folder, self.name_exp, '_rna_ribo_norm', True, self.th_out)

		logging.info('========== Get Norm Ribo-RNA : Done! ============ ')
		print ('Get Norm Ribo-RNA : Done!')

		writer.save()


	def common_to_modes(self):
		self.bam_files_info = GetPrimaryReadCountBamFiles()
		self.bam_files_info.get_all_counts(self.path_to_input_folder, self.path_to_output_folder, self.mode, self.strandedness)

		logging.info('========== Get Primary Counts : Done! ============ ')
		print ('Get Primary Counts : Done!')

		self.input_file_treatment = ReadInputFile(self.path_to_input_folder)
		self.input_file_treatment.treatment_file()
		
		logging.info('========== Treatment File : Done! ============ ')
		print ('Treatment File : Done!')

		self.reverse_translation = ReverseTranslation()
		self.reverse_translation.reverse_translation(self.input_file_treatment.peptide_mode, self.input_file_treatment.CS_mode, self.path_to_output_folder, self.name_exp)
		
		logging.info('========== Reverse Translation : Done! ============ ')
		print ('Reverse Translation : Done!')

		set_peptides = set(list(self.input_file_treatment.all_mode_peptide.keys()))

		self.alignments = Alignments(self.path_to_output_folder, self.name_exp)
		self.perfect_alignments, peptides_with_alignments = self.alignments.alignment_cs_to_genome(set_peptides)

		logging.info('========== Alignment : Done! ============ ')
		print ('Alignment : Done!')

		# positions_mcs_peptides_variants_alignment[key] = [strand, local_translation_peptide, differences_pep, info_snps, differences_ntds, [],[]]
		if len(self.input_file_treatment.manual_mode) > 0 :
			for peptide, info_peptide in self.input_file_treatment.manual_mode.items() :
				coding_sequence = info_peptide[0]
				position = info_peptide[1]
				strand = info_peptide[2]
				key = peptide+'_'+position+'_'+coding_sequence
				peptides_with_alignments.add(peptide)
				self.perfect_alignments[key] = [strand, peptide, ['NA'], ['NA'], ['NA'], [], []]

		missed_peptides = list(set_peptides - peptides_with_alignments)

		with open(self.path_to_output_folder+'alignments/missed_peptides.info', 'w') as f:
			for item in missed_peptides:
				f.write("%s\t\n" % item)

		logging.info('Total missed_peptides : %s. Find the list in : %s.', str(len(missed_peptides)), self.path_to_output_folder+'alignments/missed_peptides.info')

		logging.info('========== Common_to_modes : Done! ============ ')
		print ('common_to_modes : Done!')


	def get_annotations(self):
		info_peptide_alignments = self.get_info_peptide_alignments()

		get_info_transcripts = InfoTranscripts()
		get_info_transcripts.set_values()

		intsersect_to_annotations = IntersectAnnotations(self.perfect_alignments, self.path_to_output_folder, self.name_exp)
		intsersect_to_annotations.generate_BED_files()
		intsersect_to_annotations.perform_intersection_with_annotation()

		get_biotype = BiotypeAssignation(self.path_to_output_folder, self.name_exp, self.mode)
		
		get_biotype.get_biotypes(info_peptide_alignments, self.input_file_treatment.peptides_by_type, self.bam_files_info.bam_files_list, self.bam_files_info.bam_ribo_files_list)
		get_biotype.get_global_annotation()
		
		

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

	args = parser.parse_args()
	
	path_to_input_folder = args.path_to_input_folder
	name_exp = args.name_exp.lower()
	mode = args.mode.lower()
	strandedness = args.strandedness
	th_out = args.th_out

	if path_to_input_folder[-1] != '/':
		path_to_input_folder += '/'

	path_to_output_folder = directories_creation(path_to_input_folder, name_exp, mode, strandedness)

	t0 = time.time()

	BamQuery_obj = BamQuery(path_to_input_folder, path_to_output_folder, name_exp, mode, strandedness, th_out)
	
	t2 = time.time()
	total = t2-t0
	
	logging.info('Total time run function BamQuery to end : %f min', (total/60.0))


if __name__ == "__main__":
	main(sys.argv[1:])


#python /u/ruizma/BAM_Query/Scripts/Python/BamQuery_class.py -p /u/ruizma/BAM_Query/Res/BamQuery_exp_class/peptides.tsv -b /u/ruizma/BAM_Query/Res/BamQuery_exp_class/BAM_directories.tsv -s /u/ruizma/BAM_Query/Res/BamQuery_exp_class/ -n First_Exp
# time python /u/ruizma/BAM_Query/Scripts/Python/BamQuery_class.py -p /u/ruizma/BAM_Query/Res/BamQuery_exp_class/peptides.tsv -b /u/ruizma/BAM_Query/Res/BamQuery_exp_class/BAM_directories.tsv -s /u/ruizma/BAM_Query/Res/BamQuery_exp_class/ -n First_Exp > /u/ruizma/BAM_Query/Res/BamQuery_exp_class/test_reads.reads
# time python /u/ruizma/BAM_Query/Scripts/Python/BamQuery_class.py -p /u/ruizma/BAM_Query/Res/BamQuery_exp_class/peptides.tsv -b /u/ruizma/BAM_Query/Res/BamQuery_exp_class/BAM_directories.tsv -s /u/ruizma/BAM_Query/Res/BamQuery_exp_class/ -n First_Exp

#jupyter3
#time python /u/ruizma/BAM_Query/Scripts/Python/BamQuery/BamQuery.py -i /u/ruizma/BAM_Query/Res/BamQuery_exp_class_2/Input -n First_Exp


# time /u/ruizma/BAM_Query/BamQuery/BamQuery_normalisations_5.sh -x full -o /u/ruizma/BAM_Query/Res/BamQuery_exp/ -n test


