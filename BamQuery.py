import time, sys, os, datetime, argparse, logging, _thread, shutil, pickle
from os import listdir
from os.path import isfile, join
import pandas as pd
from pandas import ExcelFile

from readers.read_input import ReadInputFile
from readers.info_transcripts_annotation import InfoTranscripts
from readers.intersection_alignments_annotations import IntersectAnnotations

from utils.get_information_bam_files import GetInformationBamFiles
from utils.reverse_translation import ReverseTranslation
from utils.paths_arrangements import *

from genomics.alignments import *
from genomics.get_counts import GetCounts
from genomics.get_counts_sc import GetCountsSC
from genomics.normalization import Normalization
from genomics.get_biotype import BiotypeAssignation

import plotting.plots as plots

__author__ = "Maria Virginia Ruiz Cuevas"
__email__ = "maria.virginia.ruiz.cuevas@umontreal.ca"


class BamQuery:

	def __init__(self, path_to_input_folder, path_to_output_folder, name_exp, mode, strandedness, th_out, light, dev, plots, dbSNP, c, super_logger, bam_files_logger, sc, var, maxmm, genome_version):
		self.path_to_input_folder = path_to_input_folder
		self.path_to_output_folder = path_to_output_folder
		self.name_exp = name_exp
		self.strandedness = strandedness
		self.mode = mode
		self.th_out = th_out
		self.light = light
		self.dev = dev
		self.plots = plots
		self.dbSNP = dbSNP
		self.common = c
		self.super_logger = super_logger
		self.sc = sc
		self.var = var
		self.maxmm = maxmm
		self.genome_version = genome_version

		if self.mode == 'normal':
			if self.sc :
				self.run_bam_query_sc_mode(bam_files_logger)
			else:
				self.run_bam_query_normal_mode(bam_files_logger)
		else:
			self.run_bam_query_translation_mode(bam_files_logger)

		if not self.light:
			self.get_annotations()


	def run_bam_query_sc_mode(self, bam_files_logger):
		self.common_to_modes(bam_files_logger)
		
		name_path_normal = self.path_to_output_folder+'/res/'+self.name_exp+'_rna_count.csv'
		name_path_light = self.path_to_output_folder+'/res_light/'+self.name_exp+'_rna_count.csv'
		
		exists_normal = os.path.exists(name_path_normal) 
		exists_light = os.path.exists(name_path_light) 

		if not self.light:
			name_path = self.path_to_output_folder+'/res/'+self.name_exp+'_rna_count.csv'
		else:
			name_path = self.path_to_output_folder+'/res_light/'+self.name_exp+'_rna_count.csv'

		if (self.light and not exists_light) or (not self.light and not exists_normal and not exists_light):
			
			get_counts = GetCountsSC(self.path_to_output_folder, self.name_exp, self.mode, self.light, self.input_file_treatment.peptides_by_type, self.super_logger)
			
			res = get_counts.get_counts(self.perfect_alignments, self.bam_files_info.bam_files_list)
			df_counts_rna = res[0]
			self.perfect_alignments = res[1]
			df_all_alignments_rna = res[2] 

			self.super_logger.info('========== Get Norm RNA : Done! ============ ')

		elif (not self.light and not exists_normal and exists_light):
			print ('Information count already collected from light mode, filtering information for the peptides of interest !')

			self.super_logger.info('Information count already collected for light mode, filtering information for the peptides of interest !')

			name_path_light = self.path_to_output_folder+'/res_light/'+self.name_exp+'_rna_count.csv'

			df_counts_all_alignments = pd.read_csv(self.path_to_output_folder+'/res_light/'+self.name_exp+'_rna_count_All_alignments.csv', header=0, index_col=None)
			
			df_all_alignments_rna = df_counts_all_alignments[df_counts_all_alignments['Peptide'].isin(self.set_peptides) == True]
			df_all_alignments_rna = df_all_alignments_rna.set_index('Peptide')
			df_all_alignments_rna.to_csv(self.path_to_output_folder+'/res/'+self.name_exp+'_rna_count_All_alignments.csv', index=True, header=True)

			self.super_logger.info('Information All alignments for peptides of interest collected!')

			df_counts_rna_light = pd.read_csv(self.path_to_output_folder+'/res_light/'+self.name_exp+'_rna_count.csv', header=0, index_col=None)
			
			df_counts_rna = df_counts_rna_light[df_counts_rna_light['Peptides'].isin(self.set_peptides) == True]
			df_counts_rna = df_counts_rna.set_index('Peptides')
			df_counts_rna.to_csv(self.path_to_output_folder+'/res/'+self.name_exp+'_rna_count.csv', index=True, header=True)

			self.super_logger.info('Information rna counts for peptides of interest collected!')

			self.super_logger.info('Information for peptides of interest collected!')

		else:
			self.super_logger.info('Information count and normalisation already collected !')
			print ('Information count and normalisation already collected !')

	def run_bam_query_normal_mode(self, bam_files_logger):
		self.common_to_modes(bam_files_logger)
		
		name_path_normal = self.path_to_output_folder+'/res/'+self.name_exp+'_count_norm_info.xlsx'
		name_path_light = self.path_to_output_folder+'/res_light/'+self.name_exp+'_count_norm_info.xlsx'
		
		exists_normal = os.path.exists(name_path_normal) 
		exists_light = os.path.exists(name_path_light) 

		if not self.light:
			name_path = self.path_to_output_folder+'/res/'+self.name_exp+'_count_norm_info.xlsx'
		else:
			name_path = self.path_to_output_folder+'/res_light/'+self.name_exp+'_count_norm_info.xlsx'

		if (self.light and not exists_light) or (not self.light and not exists_normal and not exists_light):
			
			get_counts = GetCounts(self.path_to_output_folder, self.name_exp, self.mode, self.light, self.input_file_treatment.peptides_by_type, self.super_logger)
			
			res = get_counts.get_counts(self.perfect_alignments, self.bam_files_info.bam_files_list)
			df_counts_rna = res[0]
			self.perfect_alignments = res[1]
			df_counts_filtered = res[2] 
			order = res[3] 
			order_f = res[4] 
			df_all_alignments_rna = res[5] 

			if not self.light: 
				plots.get_heat_map(df_counts_rna, self.path_to_output_folder, self.name_exp, '_rna_counts', False, order, self.th_out)

			self.super_logger.info('========== Get Count RNA : Done! ============ ')

			normalization = Normalization(self.path_to_output_folder, self.name_exp, self.input_file_treatment.all_mode_peptide, self.mode, self.light, self.super_logger)
			def_norm_rna = normalization.get_normalization(df_counts_rna, '_rna_norm.csv')
			
			
			if not self.light: 
				plots.get_heat_map(def_norm_rna, self.path_to_output_folder, self.name_exp, '_rna_norm', True, order, self.th_out)

			writer = pd.ExcelWriter(name_path, engine='xlsxwriter')
			writer.book.use_zip64()
			df_all_alignments_rna.to_excel(writer, sheet_name='Alignments Read count RNA-seq')
			df_counts_rna.to_excel(writer, sheet_name='Read count RNA-seq by peptide')
			def_norm_rna.to_excel(writer, sheet_name='log10(RPHM) RNA-seq by peptide')
			
			writer.save()
			self.super_logger.info('========== Get Norm RNA : Done! ============ ')

		elif (not self.light and not exists_normal and exists_light):
			print ('Information count and normalisation already collected from light mode, filtering information for the peptides of interest !')

			self.super_logger.info('Information count and normalisation already collected for light mode, filtering information for the peptides of interest !')

			name_path_light = self.path_to_output_folder+'/res_light/'+self.name_exp+'_count_norm_info.xlsx'

			#df_counts_all_alignments = pd.read_csv(self.path_to_output_folder+'/res_light/temps_files/'+self.name_exp+'_rna_count_All_alignments.csv', header=0, index_col=None)
			df_counts_all_alignments = pd.read_excel(name_path_light, sheet_name='Alignments Read count RNA-seq', header=0, index_col=None, engine='openpyxl')

			df_all_alignments_rna = df_counts_all_alignments[df_counts_all_alignments['Peptide'].isin(self.set_peptides) == True]
			df_all_alignments_rna = df_all_alignments_rna.set_index('Peptide')
			df_all_alignments_rna.to_csv(self.path_to_output_folder+'/res/temps_files/'+self.name_exp+'_rna_count_All_alignments.csv', index=True, header=True)

			self.super_logger.info('Information All alignments for peptides of interest collected!')

			#df_counts_rna_light = pd.read_csv(self.path_to_output_folder+'/res_light/temps_files/'+self.name_exp+'_rna_count.csv', header=0, index_col=None)
			df_counts_rna_light = pd.read_excel(name_path_light, sheet_name='Read count RNA-seq by peptide', header=0, index_col=None, engine='openpyxl')

			df_counts_rna = df_counts_rna_light[df_counts_rna_light['Peptides'].isin(self.set_peptides) == True]
			df_counts_rna = df_counts_rna.set_index('Peptides')
			df_counts_rna.to_csv(self.path_to_output_folder+'/res/temps_files/'+self.name_exp+'_rna_count.csv', index=True, header=True)

			self.super_logger.info('Information rna counts for peptides of interest collected!')

			normalization = Normalization(self.path_to_output_folder, self.name_exp, self.input_file_treatment.all_mode_peptide, self.mode, self.light, self.super_logger)
			def_norm_rna = normalization.get_normalization(df_counts_rna, '_rna_norm.csv')
			
			self.super_logger.info('Information norm counts for peptides of interest collected!')

			writer = pd.ExcelWriter(name_path, engine='xlsxwriter')
			writer.book.use_zip64()
			df_all_alignments_rna.to_excel(writer, sheet_name='Alignments Read count RNA-seq')
			df_counts_rna.to_excel(writer, sheet_name='Read count RNA-seq by peptide')
			def_norm_rna.to_excel(writer, sheet_name='log10(RPHM) RNA-seq by peptide')
			writer.save()
			
			plots.get_heat_map(df_counts_rna, self.path_to_output_folder, self.name_exp, '_rna_counts', False, [], self.th_out)
			plots.get_heat_map(def_norm_rna, self.path_to_output_folder, self.name_exp, '_rna_norm', True, [], self.th_out)

			self.super_logger.info('Information for peptides of interest collected!')

		else:
			self.super_logger.info('Information count and normalisation already collected !')
			print ('Information count and normalisation already collected !')
		
	def run_bam_query_translation_mode(self, bam_files_logger):
		self.common_to_modes(bam_files_logger)

		if not self.light:
			name_path = self.path_to_output_folder+'/res/'+self.name_exp+'_count_norm_info.xlsx'
		else:
			name_path = self.path_to_output_folder+'/res_light/'+self.name_exp+'_count_norm_info.xlsx'
		
		exists = os.path.exists(name_path) 

		if not exists:

			get_counts = GetCounts(self.path_to_output_folder, self.name_exp, self.mode, self.light, self.input_file_treatment.peptides_by_type)
			perfect_alignments_to_return, df_counts_ribo, order, df_all_alignments_ribo = get_counts.ribo_counts(self.perfect_alignments, self.bam_files_info.bam_ribo_files_list)
			
			writer = pd.ExcelWriter(name_path, engine='xlsxwriter')
			writer.book.use_zip64()
			df_all_alignments_ribo.to_excel(writer, sheet_name='Alignments Read count Ribo-seq')
			df_counts_ribo.to_excel(writer, sheet_name='Read count Ribo-seq by peptide')

			if not self.light and len(df_counts_ribo) < 400:
				plots.get_heat_map(df_counts_ribo, self.path_to_output_folder, self.name_exp, '_ribo_counts', False, order, self.th_out)
			
			self.super_logger.info('========== Get Count Ribo : Done! ============ ')
			print ('Get Count Ribo : Done!')

			normalization = Normalization(self.path_to_output_folder, self.name_exp, self.input_file_treatment.all_mode_peptide, self.mode, self.light, self.super_logger)
			def_norm_ribo = normalization.get_normalization(df_counts_ribo, '_ribo_norm.csv')
			
			def_norm_ribo.to_excel(writer, sheet_name='log10(RPHM) Ribo-seq by peptide')
			if not self.light and len(def_norm_ribo) < 400:
				plots.get_heat_map(def_norm_ribo, self.path_to_output_folder, self.name_exp, '_ribo_norm', True, order, self.th_out)

			self.super_logger.info('========== Get Norm Ribo : Done! ============ ')
			print ('Get Norm Ribo : Done!')

			res = get_counts.get_counts(perfect_alignments_to_return, self.bam_files_info.bam_files_list)
			df_counts_rna = res[0] 
			self.perfect_alignments = res[1]
			df_counts_filtered = res[2] 
			order = res[3] 
			order_f = res[4] 
			df_all_alignments_rna = res[5] 

			df_all_alignments_rna.to_excel(writer, sheet_name='Alignments Read count RNA-seq')
			df_counts_rna.to_excel(writer, sheet_name='Read count RNA-seq by peptide')

			if not self.light and len(df_counts_rna) < 400:
				plots.get_heat_map(df_counts_rna, self.path_to_output_folder, self.name_exp, '_rna_counts', False, order, self.th_out)
			
			if not self.light and len(df_counts_filtered) < 400:
				plots.get_heat_map(df_counts_filtered, self.path_to_output_folder, self.name_exp, '_rna_ribo_counts', False, order_f, self.th_out)

			self.super_logger.info('========== Get Count RNA : Done! ============ ')
			print ('Get Count RNA : Done!')

			normalization = Normalization(self.path_to_output_folder, self.name_exp, self.input_file_treatment.all_mode_peptide, self.mode, self.light, self.super_logger)
			def_norm_rna = normalization.get_normalization(df_counts_rna, '_rna_norm.csv')
			def_norm_rna.to_excel(writer, sheet_name='log10(RPHM) RNA-seq by peptide')

			if not self.light and len(def_norm_rna) < 400:
				plots.get_heat_map(def_norm_rna, self.path_to_output_folder, self.name_exp, '_rna_norm', True, order, self.th_out)

			self.super_logger.info('========== Get Norm RNA : Done! ============ ')
			print ('Get Norm RNA : Done!')

			def_norm_rna_ribo = normalization.get_normalization(df_counts_filtered, '_rna_ribo_norm.csv')
			def_norm_rna_ribo.to_excel(writer, sheet_name='log10(RPHM) RNA and Ribo counts')
			
			if not self.light and len(def_norm_rna_ribo) < 400:
				plots.get_heat_map(def_norm_rna_ribo, self.path_to_output_folder, self.name_exp, '_rna_ribo_norm', True, order_f, self.th_out)

			self.super_logger.info('========== Get Norm Ribo-RNA : Done! ============ ')
			print ('Get Norm Ribo-RNA : Done!')

			writer.save()

		else:
			print ('Information count and normalisation already collected !')
		

	def common_to_modes(self, bam_files_logger):

		self.bam_files_info = GetInformationBamFiles(self.path_to_input_folder, self.path_to_output_folder, self.mode, self.strandedness, self.light, bam_files_logger, self.sc)

		self.super_logger.info('Total Bam Files to Query : %d.', len(self.bam_files_info.bam_files_list))

		if not self.sc :
			self.super_logger.info('========== Get Primary Counts : Done! ============ ')
			print ('Get Primary Counts : Done!')

		self.input_file_treatment = ReadInputFile(self.path_to_input_folder, self.super_logger)
		self.input_file_treatment.treatment_file()

		# if self.dev:
		# 	with open(self.path_to_output_folder+'genome_alignments/peptides_by_type.dic', 'wb') as handle:
		# 		pickle.dump(self.input_file_treatment.peptides_by_type, handle, protocol=pickle.HIGHEST_PROTOCOL)

		# 	with open(self.path_to_output_folder+'genome_alignments/all_mode_peptide.dic', 'wb') as handle:
		# 		pickle.dump(self.input_file_treatment.all_mode_peptide, handle, protocol=pickle.HIGHEST_PROTOCOL)

		# 	with open(self.path_to_output_folder+'genome_alignments/CS_mode.dic', 'wb') as handle:
		# 		pickle.dump(self.input_file_treatment.CS_mode, handle, protocol=pickle.HIGHEST_PROTOCOL)

		# 	with open(self.path_to_output_folder+'genome_alignments/peptides_mode.dic', 'wb') as handle:
		# 		pickle.dump(self.input_file_treatment.peptide_mode, handle, protocol=pickle.HIGHEST_PROTOCOL)

		# 	with open(self.path_to_output_folder+'genome_alignments/manual_mode.dic', 'wb') as handle:
		# 		pickle.dump(self.input_file_treatment.manual_mode, handle, protocol=pickle.HIGHEST_PROTOCOL)
		
		
		self.super_logger.info('========== Treatment File : Done! ============ ')
		print ('Treatment File : Done!')
		
		self.set_peptides = set(list(self.input_file_treatment.all_mode_peptide.keys()))

		if len(self.input_file_treatment.peptide_mode) > 0 or len(self.input_file_treatment.CS_mode) > 0 :

			self.reverse_translation = ReverseTranslation()
			output_message = self.reverse_translation.reverse_translation(self.input_file_treatment.peptide_mode, self.input_file_treatment.CS_mode, self.path_to_output_folder, self.name_exp)
			
			self.super_logger.info(output_message)

			self.super_logger.info('========== Reverse Translation : Done! ============ ')
			print ('Reverse Translation : Done!')
			
			#self.alignments = Alignments(self.path_to_output_folder, self.name_exp, self.light, self.dbSNP, self.common, self.super_logger, self.var, self.maxmm, self.genome_version)
			self.perfect_alignments, peptides_with_alignments = alignment_cs_to_genome(self.set_peptides, self.path_to_output_folder, self.name_exp, self.light, self.dbSNP, self.common, self.super_logger, self.var, self.maxmm, self.genome_version)

			self.super_logger.info('========== Alignment : Done! ============ ')
			print ('Alignment : Done!')

		if len(self.input_file_treatment.manual_mode) > 0 :

			not_in = False
			for peptide, info_peptide in self.input_file_treatment.manual_mode.items() :
				coding_sequence = info_peptide[0]
				position = info_peptide[1]
				strand = info_peptide[2]
				key = peptide+'_'+position+'_'+coding_sequence
				if key not in self.perfect_alignments.keys():
					not_in = True
					peptides_with_alignments.add(peptide)
					self.perfect_alignments[key] = [strand, peptide, ['NA'], ['NA'], ['NA'], [], []]
			if not_in:
				if not self.light:
					name_path = self.path_to_output_folder+'alignments/Alignments_information.dic'
				else :
					name_path = self.path_to_output_folder+'alignments/Alignments_information_light.dic'

				with open(name_path, 'wb') as handle:
					pickle.dump(self.perfect_alignments, handle, protocol=pickle.HIGHEST_PROTOCOL)

		# positions_mcs_peptides_variants_alignment[key] = [strand, local_translation_peptide, differences_pep, info_snps, differences_ntds, [],[]]
		exists = os.path.exists(self.path_to_output_folder+'alignments/missed_peptides.info')

		if not exists:
			missed_peptides = list(self.set_peptides - peptides_with_alignments)

			with open(self.path_to_output_folder+'alignments/missed_peptides.info', 'w') as f:
				for item in missed_peptides:
					f.write("%s\t\n" % item)

			self.super_logger.info('Total missed_peptides : %s. Find the list in : %s.', str(len(missed_peptides)), self.path_to_output_folder+'alignments/missed_peptides.info')

		self.super_logger.info('========== Common_to_modes : Done! ============ ')
		print ('common_to_modes : Done!')


	def get_annotations(self):

		self.super_logger.info('========== Running get_annotations ============ ')
		info_peptide_alignments = self.get_info_peptide_alignments()

		get_info_transcripts = InfoTranscripts()
		get_info_transcripts.set_values(self.genome_version)
		
		intersect_to_annotations = IntersectAnnotations(self.perfect_alignments, self.path_to_output_folder, self.name_exp, self.super_logger, self.genome_version)
		intersect_to_annotations.generate_BED_files()
		intersect_to_annotations.perform_intersection_with_annotation()

		self.super_logger.info('========== Intersect to Annotations : Done! ============ ')

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
		
		get_biotype = BiotypeAssignation(self.path_to_output_folder, self.name_exp, self.mode, list_bam_files_order_rna, list_bam_files_order_ribo, order_sample_bam_files_rna, order_sample_bam_files_ribo, self.dev, self.plots, self.super_logger, self.genome_version)
		get_biotype.get_biotypes(info_peptide_alignments, self.input_file_treatment.peptides_by_type)
		get_biotype.get_global_annotation()
		
		self.super_logger.info('========== Annotations : Done! ============ ')


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


def running_for_web(path_to_input_folder, name_exp, strandedness, genome_version, dbSNP, th_out = 8.55,):

	path_to_input_folder = path_to_input_folder
	mode = 'normal'
	strandedness = strandedness
	th_out = th_out
	light = False
	dev = True
	plots = True
	c = False
	sc = False
	var = False
	maxmm = False
	
	if dbSNP == 'dbSNP_149':
		dbSNP = 149
	elif dbSNP == 'dbSNP_151':
		dbSNP = 151
	elif dbSNP == 'dbSNP_155':
		dbSNP = 155
	else: dbSNP = 0

	if path_to_input_folder[-1] != '/':
		path_to_input_folder += '/'

	path_to_output_folder, super_logger, bam_files_logger  = directories_creation(path_to_input_folder, name_exp, mode, strandedness, light)

	t0 = time.time()

	BamQuery(path_to_input_folder, path_to_output_folder, name_exp, mode, strandedness, th_out, light, dev, plots, dbSNP, c, super_logger, bam_files_logger, sc, var, maxmm, genome_version)
	
	t2 = time.time()
	total = t2-t0
	
	# try:
	# 	super_logger.info(' ========== Removing Temporal Files ============ ')
	# 	shutil.rmtree(path_to_output_folder+'logs', ignore_errors=True)
	# 	shutil.rmtree(path_to_output_folder+'genome_alignments', ignore_errors=True)
	# 	shutil.rmtree(path_to_output_folder+'alignments', ignore_errors=True)
	# 	shutil.rmtree(path_to_output_folder+'res/BED_files', ignore_errors=True)
	# 	shutil.rmtree(path_to_output_folder+'res/AUX_files', ignore_errors=True)
	# 	shutil.rmtree(path_to_output_folder+'res/temps_files', ignore_errors=True)
	# 	os.remove(path_to_output_folder+"res/info_bam_files_tissues.csv")

	# except FileNotFoundError:
	# 	pass

	super_logger.info('Total time run function BamQuery to end : %f min', (total/60.0))
	
	return path_to_output_folder


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
	parser.add_argument('--dbSNP', type=int, default = 149,
						help='BamQuery dbSNP : 149 / 151 / 155 / 0')
	parser.add_argument('--dev', action='store_true')
	parser.add_argument('--plots', action='store_true')
	parser.add_argument('--c', action='store_true',
						help='Take into account the COMMON SNPs from the dbSNP database chosen')
	parser.add_argument('--sc', action='store_true',
						help='Query Single Cell Bam Files')
	parser.add_argument('--var', action='store_true',
						help='Keep Variants Alignments')
	parser.add_argument('--maxmm', action='store_true',
						help='Keep High Amount Alignments')
	parser.add_argument('--genome_version', type=str, default = 'v26_88',
						help='Genome version supported : v26_88 / v33_99 / v38_104')

	args = parser.parse_args()
	
	path_to_input_folder = args.path_to_input_folder
	name_exp = args.name_exp
	mode = args.mode.lower()
	dbSNP = args.dbSNP
	strandedness = args.strandedness
	th_out = args.th_out
	light = args.light
	dev = args.dev
	plots = args.plots
	c = args.c
	sc = args.sc
	var = args.var
	maxmm = args.maxmm
	genome_version = args.genome_version

	if sc :
		mode = 'normal'
		plots = False

	if (mode != 'normal' and mode != 'translation') or (dbSNP != 0 and dbSNP != 149 and dbSNP != 151 and dbSNP != 155) or (genome_version != 'v26_88' and genome_version != 'v33_99' and genome_version != 'v38_104' ):
		sys.stderr.write('error: %s\n' % 'Some arguments are not valid!')
		parser.print_help()
		sys.exit(2)

	if path_to_input_folder[-1] != '/':
		path_to_input_folder += '/'

	path_to_output_folder, super_logger, bam_files_logger = directories_creation(path_to_input_folder, name_exp, mode, strandedness, light)

	t0 = time.time()

	super_logger.info('=============== BamQuery id : %s ===============', name_exp)
	super_logger.info('=============== Parameters ===============')
	super_logger.info(' - Mode : %s , Strandedness :  %s, Light:  %s ', mode, strandedness, str(light) )
	super_logger.info(' - Single-Cell experiment (sc) :  %s', str(sc))
	super_logger.info(' - dbSNP :  %s, COMMON SNPs : %s, Genome Version : %s ', str(dbSNP), str(c), genome_version)
	super_logger.info(' - Plots : %s', str(plots))
	super_logger.info(' - Keep Variant Alignments : %s, Keep High Amount Alignments : %s', str(var), str(maxmm))

	BamQuery(path_to_input_folder, path_to_output_folder, name_exp, mode, strandedness, th_out, light, dev, plots, dbSNP, c, super_logger, bam_files_logger, sc, var, maxmm, genome_version)
	

	if not dev:
		
		if not light:
			try:
				shutil.rmtree(path_to_output_folder+'genome_alignments')
				shutil.rmtree(path_to_output_folder+'res/BED_files')
				shutil.rmtree(path_to_output_folder+'res/AUX_files')
				shutil.rmtree(path_to_output_folder+'res/temps_files')
			except FileNotFoundError:
				pass
		else:
			try:
				shutil.rmtree(path_to_output_folder+'res_light/AUX_files')
				shutil.rmtree(path_to_output_folder+'res_light/temps_files')
				shutil.rmtree(path_to_output_folder+'res_light/BED_files')
			except FileNotFoundError:
				pass

	t2 = time.time()
	total = t2-t0
	
	super_logger.info('Total time run function BamQuery to end : %f min', (total/60.0))


if __name__ == "__main__":
	main(sys.argv[1:])


