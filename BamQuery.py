import time, sys, os, argparse, logging, shutil, pickle
import pandas as pd

from readers.read_input import ReadInputFile
from readers.intersection_alignments_annotations import IntersectAnnotations

from utils.get_information_bam_files import GetInformationBamFiles
from utils.reverse_translation import ReverseTranslation
from utils.immunogenicity import Immunogenicity
from utils.paths_arrangements import *

from genomics.alignments import *
from genomics.get_counts import GetCounts
from genomics.get_counts_sc import GetCountsSC
from genomics.normalization import Normalization
from genomics.get_biotype import BiotypeAssignation

import plotting.plots as plots

__author__ = "Maria Virginia Ruiz Cuevas"
__email__ = "maria.virginia.ruiz.cuevas@umontreal.ca"

path_to_lib = '/'.join(os.path.abspath(__file__).split('/')[:-2])+'/lib/'

class BamQuery:

	def __init__(self, path_to_input_folder, path_to_output_folder, name_exp, mode, strandedness, th_out, light, dev, plots, dbSNP, c, super_logger, bam_files_logger, sc, var, maxmm, genome_version, overlap, mouse, t):
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
		self.overlap =  overlap
		self.mouse = mouse
		self.threads = t

		if self.mode == 'normal':
			if self.sc :
				self.run_bam_query_sc_mode(bam_files_logger)
			else:
				self.run_bam_query_normal_mode(bam_files_logger)
				if not self.light:
					self.get_annotations()
		else:
			self.run_bam_query_translation_mode(bam_files_logger)



	def run_bam_query_sc_mode(self, bam_files_logger):
		self.common_to_modes(bam_files_logger)
		
		name_path_normal = self.path_to_output_folder+'res/'+self.name_exp+'_rna_count.csv'
		
		exists_normal = os.path.exists(name_path_normal) 
		if not exists_normal :
			
			get_counts = GetCountsSC(self.path_to_output_folder, self.name_exp, self.mode, self.light, self.input_file_treatment.all_mode_peptide, self.super_logger, self.threads)
			res = get_counts.get_counts(self.perfect_alignments, self.bam_files_info.bam_files_list)
			df_counts_rna = res[0]
			self.perfect_alignments = res[1]
			df_all_alignments_rna = res[2] 

			self.super_logger.info('========== Get Count RNA single cell : Done! ============ ')


	def run_bam_query_normal_mode(self, bam_files_logger):
		self.common_to_modes(bam_files_logger)
		
		name_path_normal = self.path_to_output_folder+'res/'+self.name_exp+'_count_norm_info.xlsx'
		name_path_light = self.path_to_output_folder+'res_light/'+self.name_exp+'_count_norm_info.xlsx'
		
		exists_normal = os.path.exists(name_path_normal) 
		exists_light = os.path.exists(name_path_light) 

		if not self.light:
			name_path = self.path_to_output_folder+'res/'+self.name_exp+'_count_norm_info.xlsx'
			path_temps_file = self.path_to_output_folder+'res/temps_files/'
		else:
			name_path = self.path_to_output_folder+'res_light/'+self.name_exp+'_count_norm_info.xlsx'
			path_temps_file = self.path_to_output_folder+'res_light/temps_files/'

		if (self.light and not exists_light) or (not self.light and not exists_normal and not exists_light):
			
			get_counts = GetCounts(self.path_to_output_folder, self.name_exp, self.mode, self.light,  self.input_file_treatment.all_mode_peptide, self.super_logger, self.threads)
			res = get_counts.get_counts(self.perfect_alignments, self.bam_files_info.bam_files_list, self.overlap)
			df_counts_rna = res[0]
			self.perfect_alignments = res[1]
			df_all_alignments_rna = res[2] 

			if not self.light: 
				plots.get_heat_map(df_counts_rna, self.path_to_output_folder+'plots/heat_maps/transcription_evidence_heatmap/', self.mode, path_temps_file, self.name_exp, '_rna_counts', False, self.th_out)

			self.super_logger.info('========== Get Count RNA : Done! ============ ')

			normalization = Normalization(self.path_to_output_folder, self.name_exp, self.input_file_treatment.all_mode_peptide, self.mode, self.light, self.super_logger, self.dev)
			def_norm_rna = normalization.get_normalization(df_counts_rna, '_rna_norm.csv')
			
			if not self.light: 
				plots.get_heat_map(def_norm_rna, self.path_to_output_folder+'plots/heat_maps/transcription_evidence_heatmap/', self.mode, path_temps_file, self.name_exp, '_rna_norm', True, self.th_out)

			df_counts_rna.reset_index(inplace=True)
			writer = pd.ExcelWriter(name_path, engine='xlsxwriter')
			writer.book.use_zip64()

			if len(df_all_alignments_rna) < 1048576 and os.path.getsize(path_temps_file+self.name_exp+'_rna_count_All_alignments.csv') < 700000000: 
				df_all_alignments_rna.to_excel(writer, sheet_name='Alignments Read count RNA-seq',index=False)
			else:
				if self.light:
					df_all_alignments_rna.to_csv(self.path_to_output_folder+'res_light/'+self.name_exp+'_rna_count_All_alignments.csv', index=False)
				else:
					df_all_alignments_rna.to_csv(self.path_to_output_folder+'res/'+self.name_exp+'_rna_count_All_alignments.csv', index=False)
				
			df_counts_rna.to_excel(writer, sheet_name='Read count RNA-seq by peptide',index=False)
			def_norm_rna.to_excel(writer, sheet_name='log10(RPHM) RNA-seq by peptide',index=False)
			
			writer.save()
			self.super_logger.info('========== Get Norm RNA : Done! ============ ')

		elif (not self.light and not exists_normal and exists_light):
			print ('Information count and normalisation already collected from light mode, filtering information for the peptides of interest !')

			self.super_logger.info('Information count and normalisation already collected for light mode, filtering information for the peptides of interest !')

			name_path_light = self.path_to_output_folder+'/res_light/'+self.name_exp+'_count_norm_info.xlsx'

			df_counts_all_alignments = pd.read_excel(name_path_light, sheet_name='Alignments Read count RNA-seq', header=0, index_col=False, engine='openpyxl')

			df_all_alignments_rna =  [df_counts_all_alignments['Peptide'].isin(self.set_peptides)]
			df_all_alignments_rna.to_csv(self.path_to_output_folder+'/res/temps_files/'+self.name_exp+'_rna_count_All_alignments.csv', index=False, header=True)

			self.super_logger.info('Information All alignments for peptides of interest collected!')

			df_counts_rna_light = pd.read_excel(name_path_light, sheet_name='Read count RNA-seq by peptide', header=0, index_col=False, engine='openpyxl')
			df_counts_rna = df_counts_rna_light[df_counts_rna_light['Peptide'].isin(self.set_peptides)]
			
			df_counts_rna.to_csv(self.path_to_output_folder+'/res/temps_files/'+self.name_exp+'_rna_count.csv', index=False, header=True)

			self.super_logger.info('Information rna counts for peptides of interest collected!')

			normalization = Normalization(self.path_to_output_folder, self.name_exp, self.input_file_treatment.all_mode_peptide, self.mode, self.light, self.super_logger, self.dev)
			def_norm_rna = normalization.get_normalization(df_counts_rna, '_rna_norm.csv')
			
			self.super_logger.info('Information norm counts for peptides of interest collected!')

			df_counts_rna.reset_index(inplace=True)
			writer = pd.ExcelWriter(name_path, engine='xlsxwriter')
			writer.book.use_zip64()
			if len(df_all_alignments_rna) < 1048576:
				df_all_alignments_rna.to_excel(writer, sheet_name='Alignments Read count RNA-seq',index=False)
			else:
				df_all_alignments_rna.to_csv(writer, index=False, header = 0)
			df_counts_rna.to_excel(writer, sheet_name='Read count RNA-seq by peptide',index=False)
			def_norm_rna.to_excel(writer, sheet_name='log10(RPHM) RNA-seq by peptide',index=False)
			writer.save()
			
			plots.get_heat_map(df_counts_rna, self.path_to_output_folder+'plots/heat_maps/transcription_evidence_heatmap/', self.mode, path_temps_file, self.name_exp, '_rna_counts', False, self.th_out)
			plots.get_heat_map(def_norm_rna, self.path_to_output_folder+'plots/heat_maps/transcription_evidence_heatmap/', self.mode, path_temps_file, self.name_exp, '_rna_norm', True, self.th_out)

			self.super_logger.info('Information for peptides of interest collected!')

		else:
			self.super_logger.info('Information count and normalisation already collected !')
			print ('Information count and normalisation already collected !')
			
			if self.light:
				name_path = self.path_to_output_folder +'/alignments/Alignments_information_light_rna.dic'
			else:
				name_path = self.path_to_output_folder +'/alignments/Alignments_information_rna.dic'
				exists_alignments_information_rna = os.path.exists(name_path) 
				if not exists_alignments_information_rna:
					name_path = self.path_to_output_folder +'/alignments/Alignments_information_light_rna.dic'
			
			with open(name_path, 'rb') as fp:
				self.perfect_alignments = pickle.load(fp)
		
	def run_bam_query_translation_mode(self, bam_files_logger):
		self.common_to_modes(bam_files_logger)

		path_temps_file = self.path_to_output_folder+'/res_translation/temps_files'
		name_path = self.path_to_output_folder+'/res_translation/'+self.name_exp+'_ribo_count_info.xlsx'
		exists = os.path.exists(name_path) 

		if not exists:
			get_counts = GetCounts(self.path_to_output_folder, self.name_exp, self.mode, self.light, self.input_file_treatment.all_mode_peptide, self.super_logger, self.threads)
			res = get_counts.get_counts(self.perfect_alignments, self.bam_files_info.bam_files_list, True)
			df_counts_ribo = res[0]
			self.perfect_alignments = res[1]
			df_all_alignments_ribo = res[2] 

			if not self.light: 
				plots.get_heat_map(df_counts_ribo, self.path_to_output_folder+'plots/heat_maps/translation_evidence_heatmap/', self.mode, path_temps_file, self.name_exp, '_ribo_counts', False, self.th_out)

			self.super_logger.info('========== Get Count Ribo : Done! ============ ')

			normalization = Normalization(self.path_to_output_folder, self.name_exp, self.input_file_treatment.all_mode_peptide, self.mode, self.light, self.super_logger, self.dev)
			def_norm_ribo = normalization.get_normalization(df_counts_ribo, '_ribo_norm.csv')
			
			if not self.light: 
				plots.get_heat_map(def_norm_ribo, self.path_to_output_folder+'plots/heat_maps/translation_evidence_heatmap/', self.mode, path_temps_file, self.name_exp, '_ribo_norm', True, self.th_out)

			df_counts_ribo.reset_index(inplace=True)
			writer = pd.ExcelWriter(name_path, engine='xlsxwriter')
			writer.book.use_zip64()
			df_all_alignments_ribo.to_excel(writer, sheet_name='Alignments Read count Ribo-seq',index=False)
			df_counts_ribo.to_excel(writer, sheet_name='Read count Ribo-seq by peptide',index=False)
			def_norm_ribo.to_excel(writer, sheet_name='log10(RPHM) Ribo-seq by peptide',index=False)
			
			writer.save()
			self.super_logger.info('========== Get Norm Ribo : Done! ============ ')
			

	def common_to_modes(self, bam_files_logger):

		self.bam_files_info = GetInformationBamFiles(self.path_to_input_folder, self.path_to_output_folder, self.mode, self.strandedness, self.light, bam_files_logger, self.sc, self.genome_version, self.mouse, self.threads)

		bam_files_logger 
		handlers = bam_files_logger.handlers[:]
		for handler in handlers:
			bam_files_logger.removeHandler(handler)
			handler.close()

		self.super_logger.info('Total Bam Files to Query : %d.', len(self.bam_files_info.bam_files_list))

		self.input_file_treatment = ReadInputFile(self.path_to_input_folder, self.super_logger, self.genome_version)
		self.input_file_treatment.treatment_file()

		if self.dev:
		 	with open(self.path_to_output_folder+'genome_alignments/peptides_by_type_user.dic', 'wb') as handle:
		 		pickle.dump(self.input_file_treatment.peptides_by_type_user, handle, protocol=pickle.HIGHEST_PROTOCOL)

		self.super_logger.info('========== Treatment File : Done! ============ ')
		print ('Treatment File : Done!')
		
		self.set_peptides = set(list(self.input_file_treatment.all_mode_peptide.keys()))

		if self.mouse:
			if self.genome_version == 'M24':
				index_genome = path_to_lib+'genome_versions/genome_mouse_m24/Index_STAR_2.7.9a/'
				chrs_info = index_genome+'/chrName.txt'
				splice_junctions  = index_genome + 'sjdbList.fromGTF.out.tab'

			if self.genome_version == 'M30':
				index_genome = path_to_lib+'genome_versions/genome_mouse_m30/Index_STAR_2.7.9a/'
				chrs_info = index_genome+'/chrName.txt'
				splice_junctions  = index_genome + 'sjdbList.fromGTF.out.tab'
		else:
			if self.genome_version == 'v26_88': 
				index_genome = path_to_lib+'genome_versions/genome_v26_88/Index_STAR_2.7.9a/'
				chrs_info = index_genome+'/chrName.txt'
				splice_junctions  = index_genome + 'sjdbList.fromGTF.out.tab'

			elif self.genome_version == 'v33_99':
				index_genome = path_to_lib+'genome_versions/genome_v33_99/Index_STAR_2.7.9a/'
				chrs_info = index_genome+'/chrName.txt'
				splice_junctions  = index_genome + 'sjdbList.fromGTF.out.tab'

			else:
				index_genome = path_to_lib+'genome_versions/genome_v38_104/Index_STAR_2.7.9a/'
				chrs_info = index_genome+'/chrName.txt'
				splice_junctions  = index_genome + 'sjdbList.fromGTF.out.tab'

		chrs_info = pd.read_csv(chrs_info, header=None)
		references = list(chrs_info[0])
		path_to_save_list = self.path_to_output_folder+'/genome_alignments/references_chrs.pkl'
		
		with open(path_to_save_list, "wb") as f:
			pickle.dump(references, f)

		if len(self.input_file_treatment.peptide_mode) > 0 or len(self.input_file_treatment.CS_mode) > 0 :

			self.reverse_translation = ReverseTranslation()
			output_message = self.reverse_translation.reverse_translation(self.input_file_treatment.peptide_mode, self.input_file_treatment.CS_mode, self.path_to_output_folder, self.name_exp, self.threads)
			
			self.super_logger.info(output_message)

			self.super_logger.info('========== Reverse Translation : Done! ============ ')
			print ('Reverse Translation : Done!')
			
			#self.alignments = Alignments(self.path_to_output_folder, self.name_exp, self.light, self.dbSNP, self.common, self.super_logger, self.var, self.maxmm, self.genome_version)
			self.perfect_alignments, peptides_with_alignments = alignment_cs_to_genome(self.set_peptides, self.path_to_output_folder, self.name_exp, self.light, self.dbSNP, self.common, self.super_logger, self.var, self.maxmm, self.genome_version, self.mode, self.mouse, self.threads)

			self.super_logger.info('========== Alignment : Done! ============ ')
			print ('Alignment : Done!')
			if len(self.perfect_alignments) == 0:
				self.super_logger.info('========== No genomic locations were found for the peptides queried. ============')
				self.super_logger.info('========== BamQuery status: finished! ============')
				sys.exit(2)

		if len(self.input_file_treatment.manual_mode) > 0 :

			splice_junctions_annotated = pd.read_csv(splice_junctions, header=None, sep='\t')
			try:
				perfect_alignments_exists = isinstance(self.perfect_alignments, dict)
			except AttributeError:
				perfect_alignments_exists = False
				self.perfect_alignments = {}
				peptides_with_alignments = set()
			info_to_add = []

			for peptide, info_peptide in self.input_file_treatment.manual_mode.items() :
				for info in info_peptide:
					coding_sequence = info[0]
					position = info[1]
					strand = info[2]
					key = peptide+'_'+position+'_'+coding_sequence
					known_splice_junction = []
					chr = position.split(':')[0]

					if '|' in position:
						result = re.findall(r"\d+", position.split(':')[1])[1:-1]
						tuples = [(int(result[i]), int(result[i+1])) for i in range(0, len(result), 2)]
							
						for tuple in tuples:
							annotated_sj = splice_junctions_annotated[(splice_junctions_annotated[0]==chr) & (splice_junctions_annotated[1]==tuple[0]+1)& (splice_junctions_annotated[2]==tuple[1]-1) & (splice_junctions_annotated[3]==strand)]
							if not annotated_sj.empty:
								known_splice_junction.append('yes')
							else:
								known_splice_junction.append('no')
					else:
						known_splice_junction.append('NA')

					known_splice_junction = '/'.join(known_splice_junction)
					
					if key not in self.perfect_alignments.keys():
						peptides_with_alignments.add(peptide)
						self.perfect_alignments[key] = [strand, peptide, ['NA'], ['NA'], ['NA'], []]
					
					info_to_add.append([peptide, strand, position, known_splice_junction, coding_sequence, peptide])

			if not perfect_alignments_exists:
				if not self.light:
					name_path = self.path_to_output_folder+'alignments/Alignments_information.dic'
				else :
					name_path = self.path_to_output_folder+'alignments/Alignments_information_light.dic'

				with open(name_path, 'wb') as handle:
					pickle.dump(self.perfect_alignments, handle, protocol=pickle.HIGHEST_PROTOCOL)

			columns = ["Peptide", "Strand", "Alignment", 'Known Splice Junction', "MCS", "Peptide in Reference"]
			path = self.path_to_output_folder+'alignments/alignments_summary_information.pkl'
			try:
				df1 = pd.read_pickle(path)
				df_aux = pd.DataFrame(info_to_add, columns=columns)
				alignments_summary_information = pd.concat([df1, df_aux])
			except FileNotFoundError:
				alignments_summary_information = pd.DataFrame(info_to_add, columns=columns)
				
			alignments_summary_information.to_pickle(path)

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

		intersect_to_annotations = IntersectAnnotations(self.perfect_alignments, self.path_to_output_folder, self.mode, self.name_exp, self.super_logger, self.genome_version, self.mouse)
		intersect_to_annotations.generate_BED_files()
		if not self.mouse:
			intersect_to_annotations.perform_intersection_with_annotation()
		else:
			intersect_to_annotations.perform_intersection_with_annotation_mouse()

		self.super_logger.info('========== Intersect to Annotations : Done! ============ ')

		list_bam_files_order_rna = []
		order_sample_bam_files_rna = {}

		for name_sample, info_bam in sorted(self.bam_files_info.bam_files_list.items(), key=lambda e: e[1][-2], reverse=False):
			list_bam_files_order_rna.append(name_sample)
			group = info_bam[3]
			try:
				order_sample_bam_files_rna[group].append(name_sample)
			except KeyError:
				order_sample_bam_files_rna[group] = [name_sample]

		get_biotype = BiotypeAssignation(self.path_to_output_folder, self.name_exp, self.mode, list_bam_files_order_rna, order_sample_bam_files_rna, self.dev, self.plots, self.super_logger, self.genome_version, self.mouse, self.threads)
		get_biotype.get_biotypes(info_peptide_alignments, self.input_file_treatment.peptides_by_type_user)
		
		try:
			get_biotype.get_global_annotation()
			self.super_logger.info('========== Annotations : Done! ============ ')
		except MemoryError:
			self.super_logger.info('Biotype classification has stopped due to lack of memory. Please try again by allocating more memory to the process or processing fewer peptides.')
			return

	def get_info_peptide_alignments(self):

		info_peptide_alignments = {}

		for peptide_alignment in self.perfect_alignments:
			peptide = peptide_alignment.split('_')[0]
			alignment = peptide_alignment.split('_')[1]
			MCS = peptide_alignment.split('_')[2]
			count_rna = self.perfect_alignments[peptide_alignment][-1]
			strand = self.perfect_alignments[peptide_alignment][0]
			
			try:
				info_peptide = info_peptide_alignments[peptide]
				info_peptide[0][peptide_alignment] = [strand, count_rna]
				info_peptide[1] += sum(count_rna) 
			except KeyError:
				dic = {}
				dic[peptide_alignment] = [strand, count_rna]
				info_peptide_alignments[peptide] = [dic, sum(count_rna)]

		return info_peptide_alignments


def running_for_web(path_to_input_folder, name_exp, strandedness, genome_version, dbSNP, th_out = 8.55,):

	path_to_input_folder = path_to_input_folder
	mode = 'normal'
	strandedness = strandedness
	th_out = th_out
	light = False
	dev = False
	plots = True
	c = False
	sc = False
	var = False
	maxmm = False
	overlap = False
	mouse = False
	threads = 8
	
	if dbSNP == 'dbSNP_149':
		dbSNP = 149
	elif dbSNP == 'dbSNP_151':
		dbSNP = 151
	elif dbSNP == 'dbSNP_155':
		dbSNP = 155
	else: dbSNP = 0

	if path_to_input_folder[-1] != '/':
		path_to_input_folder += '/'

	path_to_output_folder, super_logger, bam_files_logger, handler_super_logger, handler_bam_files_logger  = directories_creation(path_to_input_folder, name_exp, mode, light, sc)

	t0 = time.time()

	BamQuery(path_to_input_folder, path_to_output_folder, name_exp, mode, strandedness, th_out, light, dev, plots, dbSNP, c, super_logger, bam_files_logger, sc, var, maxmm, genome_version, overlap, mouse, threads)
	
	predictions = Immunogenicity(path_to_output_folder, name_exp)
	predictions.get_predictions()

	t2 = time.time()
	total = t2-t0
	
	super_logger.info('Total time run function BamQuery to end : %f min', (total/60.0))
	handler_super_logger.close()
	handler_bam_files_logger.close()
	logging.shutdown()
	del super_logger
	del bam_files_logger
	del handler_super_logger
	del handler_bam_files_logger
	
	try:
		os.remove(path_to_output_folder+"res/info_bam_files_tissues.csv")
		shutil.rmtree(path_to_output_folder+'genome_alignments', ignore_errors=True)
		shutil.rmtree(path_to_output_folder+'alignments', ignore_errors=True)
		shutil.rmtree(path_to_output_folder+'res/BED_files', ignore_errors=True)
		shutil.rmtree(path_to_output_folder+'res/AUX_files', ignore_errors=True)
		shutil.rmtree(path_to_output_folder+'res/temps_files', ignore_errors=True)
		os.system('rmdir /S /Q "{}"'.format(path_to_output_folder+'logs'))
		os.system('rm -rf "{}"'.format(path_to_output_folder+'logs'))
		shutil.rmtree(path_to_output_folder+'logs', ignore_errors=True)	
	except Exception as err:
		print (err)
		pass
	
	path_to_readme_file = path_to_lib+'README.txt'
	shutil.copy2(path_to_readme_file, path_to_output_folder)

	return path_to_output_folder


def main(argv):

	parser = argparse.ArgumentParser(description='======== BamQuery ========')
	parser.add_argument('path_to_input_folder', type=str, help='Path to the input folder where to find BAM_directories.tsv and peptides.tsv')
	parser.add_argument('name_exp', type=str, help='BamQuery search Id')
	parser.add_argument('genome_version', type=str, help='Genome human releases : v26_88 / v33_99 / v38_104; Genome mouse releases : M24 / M30')
	parser.add_argument('--mode', type=str, default = 'normal', help='BamQuery search mode : normal / translation')
	parser.add_argument('--th_out', type=float, default = 8.55, help='Threshold to assess expression comparation with other tissues')
	parser.add_argument('--dbSNP', type=int, default = 0, help='Human dbSNP : 149 / 151 / 155')
	parser.add_argument('--c', action='store_true', help='Take into account the only common SNPs from the dbSNP database chosen')
	parser.add_argument('--strandedness', action='store_true', help='Take into account strandedness of the samples')
	parser.add_argument('--light', action='store_true', help='Display only the count and norm count for peptides and regions')
	parser.add_argument('--sc', action='store_true', help='Query Single Cell Bam Files')
	parser.add_argument('--var', action='store_true', help='Keep Variants Alignments')
	parser.add_argument('--maxmm', action='store_true', help='Enable STAR to generate a larger number of alignments')
	parser.add_argument('--overlap', action='store_true', help='Count overlapping reads')
	parser.add_argument('--plots', action='store_true', help='Plot biotype pie-charts')
	parser.add_argument('--m', action='store_true', help='Mouse genome')
	parser.add_argument('--dev', action='store_true', help='Save all temps files')
	parser.add_argument('--t', type=int, default = 4, help='Specify the number of processing threads to run BamQuery. The default is 4')
	
	args = parser.parse_args()
	path_to_input_folder = args.path_to_input_folder
	name_exp = args.name_exp
	genome_version = args.genome_version
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
	overlap = args.overlap
	mouse = args.m
	t = args.t

	if sc and mouse:
		sys.stderr.write('error: %s\n' % 'Some arguments are not valid! Please verify the use of a single BamQuery method to perform the search. (sc or mouse)')
		parser.print_help()
		sys.exit(2)

	if light:
		plots = False

	if mouse :
		c = False

	if sc :
		mode = 'normal'
		plots = False

	if (mode != 'normal' and mode != 'translation') or (dbSNP != 0 and dbSNP != 149 and dbSNP != 151 and dbSNP != 155) or (genome_version != 'v26_88' and genome_version != 'v33_99' and genome_version != 'v38_104' and genome_version != 'M24' and genome_version != 'M30' ):
		sys.stderr.write('error: %s\n' % 'Some arguments are not valid!')
		parser.print_help()
		sys.exit(2)

	if path_to_input_folder[-1] != '/':
		path_to_input_folder += '/'

	path_to_output_folder, super_logger, bam_files_logger, handler_super_logger, handler_bam_files_logger = directories_creation(path_to_input_folder, name_exp, mode, light, sc)

	t0 = time.time()

	if mouse:
		if (genome_version != 'M24' and genome_version != 'M30') or (genome_version == 'M24'):
			genome_version = 'M24'
			dbSNP = 'mouse_GRCm38'
		if genome_version == 'M30':
			genome_version = 'M30'
			dbSNP = 'mouse_GRCm39'
		super_logger.info('=============== Start Parameters ===============')
		super_logger.info(' - BamQuery id : %s ', name_exp)
		super_logger.info(' - Mode : %s, Strandedness :  %s, Light:  %s ', mode, strandedness, str(light) )
		super_logger.info(' - Single-Cell experiment (sc) :  %s', str(sc))
		super_logger.info(' - dbSNP :  %s, COMMON SNPs : %s, Genome Version : %s ', str(dbSNP), str(c), genome_version)
		super_logger.info(' - Plots : %s', str(plots))
		super_logger.info(' - Keep Variant Alignments : %s, Keep High Amount Alignments : %s', str(var), str(maxmm))
		super_logger.info(' - Counting overlapping reads : %s', str(overlap))
		super_logger.info(' - Mouse Genome : %s', str(mouse))
		super_logger.info(' - Threads : %s', str(t))
		super_logger.info('=============== End Parameters ===============')

	else:
		super_logger.info('=============== Start Parameters ===============')
		super_logger.info(' - BamQuery id : %s ', name_exp)
		super_logger.info(' - Mode : %s, Strandedness :  %s, Light:  %s ', mode, strandedness, str(light) )
		super_logger.info(' - Single-Cell experiment (sc) :  %s', str(sc))
		super_logger.info(' - dbSNP :  %s, COMMON SNPs : %s, Genome Version : %s ', str(dbSNP), str(c), genome_version)
		super_logger.info(' - Plots : %s', str(plots))
		super_logger.info(' - Keep Variant Alignments : %s, Keep High Amount Alignments : %s', str(var), str(maxmm))
		super_logger.info(' - Counting overlapping reads : %s', str(overlap))
		super_logger.info(' - Mouse Genome : %s', str(mouse))
		super_logger.info(' - Threads : %s', str(t))
		super_logger.info('=============== End Parameters ===============')

	BamQuery(path_to_input_folder, path_to_output_folder, name_exp, mode, strandedness, th_out, light, dev, plots, dbSNP, c, super_logger, bam_files_logger, sc, var, maxmm, genome_version, overlap, mouse, t)
	
	print ('========== BamQuery : Done! ============ ')
	if not dev:
		try:
			os.remove(path_to_output_folder+"alignments/Alignments_information.dic")
			os.remove(path_to_output_folder+"alignments/info_treated_bam_files.pkl")
		except:
			pass
		os.remove(path_to_output_folder+"alignments/alignments_summary_information.pkl")
		shutil.rmtree(path_to_output_folder+'genome_alignments')
		
		if sc:
			shutil.rmtree(path_to_output_folder+'res/temps_files')
			os.remove(path_to_output_folder+"alignments/Alignments_information_sc.dic")
			
		if mode == 'translation':
			shutil.rmtree(path_to_output_folder+'res_translation/BED_files')
			shutil.rmtree(path_to_output_folder+'res_translation/temps_files')
			shutil.rmtree(path_to_output_folder+'res_translation/AUX_files')
			os.remove(path_to_output_folder+"alignments/Alignments_information_ribo.dic")
		
		if mode == 'normal' and not light and not sc:
			shutil.rmtree(path_to_output_folder+'res/BED_files')
			shutil.rmtree(path_to_output_folder+'res/AUX_files')
			shutil.rmtree(path_to_output_folder+'res/temps_files')
			os.remove(path_to_output_folder+"alignments/Alignments_information_rna.dic")
		
		if light :
			shutil.rmtree(path_to_output_folder+'res_light/temps_files')
			shutil.rmtree(path_to_output_folder+'res_light/AUX_files')
			os.remove(path_to_output_folder+"alignments/Alignments_information_light_rna.dic")
			os.remove(path_to_output_folder+"alignments/Alignments_information_light.dic")


	t2 = time.time()
	total = t2-t0
	
	super_logger.info('Total time run function BamQuery to end : %f min', (total/60.0))
	logging.shutdown()

if __name__ == "__main__":
	main(sys.argv[1:])


