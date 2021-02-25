import warnings, time, sys, os, datetime, getopt, logging, _thread
warnings.filterwarnings("ignore")
from os import listdir
from os.path import isfile, join
import pandas as pd

from readers.read_input import ReadInputFile
from readers.info_transcripts_annotation import InfoTranscripts
from readers.intersection_alignments_annotations import IntersectAnnotations

from utils.primary_read_count import GetPrimaryReadCountBamFiles
from utils.reverse_translation import ReverseTranslation

from genomics.alignments import Alignments
from genomics.get_counts import GetCounts
from genomics.normalization import Normalization
from genomics.get_biotype import BiotypeAssignation

import plotting.plots as plots

__author__ = "Maria Virginia Ruiz Cuevas"


class BamQuery:

	def __init__(self, path_to_input_folder, path_to_output_folder, name_exp, mode):
		self.path_to_input_folder = path_to_input_folder
		self.path_to_output_folder = path_to_output_folder
		self.name_exp = name_exp
		self.mode = mode

		if self.mode == 'normal':
			self.run_bam_query_normal_mode()
		else:
			self.run_bam_query_filter_mode()

		self.get_annotations()

	def run_bam_query_normal_mode(self):
		self.common_to_modes()
		
		get_counts = GetCounts(self.path_to_output_folder, self.name_exp, self.mode)
		df_counts, self.perfect_alignments = get_counts.get_counts(self.res_star, self.bam_files_info.bam_files_list)
		plots.get_heat_map(df_counts, self.path_to_output_folder, self.name_exp, 'get_counts')

		normalization = Normalization(self.path_to_output_folder, self.name_exp, self.bam_files_info.bam_files_list, self.mode)
		def_norm = normalization.get_normalization(df_counts, '_norm.csv')
		plots.get_heat_map(def_norm, self.path_to_output_folder, self.name_exp, 'get_norm')


	def run_bam_query_filter_mode(self):
		self.common_to_modes()

		get_counts = GetCounts(self.path_to_output_folder, self.name_exp, self.mode)
		perfect_alignments_to_return, df_counts = get_counts.filter_counts(self.res_star, self.bam_files_info.bam_ribo_files_list)
		plots.get_heat_map(df_counts, self.path_to_output_folder, self.name_exp, 'filter_counts')
		
		normalization = Normalization(self.path_to_output_folder, self.name_exp, self.bam_files_info.bam_ribo_files_list, self.mode)
		def_norm = normalization.get_normalization(df_counts, '_filter_norm.csv')
		plots.get_heat_map(def_norm, self.path_to_output_folder, self.name_exp, 'filter_norm')

		df_counts, self.perfect_alignments = get_counts.get_counts(perfect_alignments_to_return, self.bam_files_info.bam_files_list)
		plots.get_heat_map(df_counts, self.path_to_output_folder, self.name_exp, '_counts')

		normalization = Normalization(self.path_to_output_folder, self.name_exp, self.bam_files_info.bam_files_list, self.mode)
		def_norm = normalization.get_normalization(df_counts, '_norm.csv')
		plots.get_heat_map(def_norm, self.path_to_output_folder, self.name_exp, '_norm')


	def common_to_modes(self):
		self.bam_files_info = GetPrimaryReadCountBamFiles()
		self.bam_files_info.get_all_counts(self.path_to_input_folder, self.path_to_output_folder)

		self.input_file_treatment = ReadInputFile(self.path_to_input_folder)
		self.input_file_treatment.treatment_file()
		
		self.reverse_translation = ReverseTranslation()
		self.reverse_translation.reverse_translation(self.input_file_treatment.peptide_mode, self.input_file_treatment.CS_mode, self.path_to_output_folder, self.name_exp)
		
		self.alignments = Alignments(self.path_to_output_folder, self.name_exp)
		self.res_star = self.alignments.alignment_cs_to_genome()

		if len(self.input_file_treatment.manual_mode) > 0 :
			for peptide, info_peptide in self.input_file_treatment.manual_mode.items() :
				coding_sequence = info_peptide[0]
				position = info_peptide[1]
				strand = info_peptide[2]
				key = peptide+'_'+position
				self.res_star[key] = [strand, coding_sequence, peptide, ['Peptide Manual Mode'],0]


	def get_annotations(self):
		get_info_transcripts = InfoTranscripts()
		get_info_transcripts.set_values()

		intsersect_to_annotations = IntersectAnnotations(self.perfect_alignments, self.path_to_output_folder, self.name_exp)
		intsersect_to_annotations.generate_BED_files()
		intsersect_to_annotations.perform_intersection_with_annotation()

		get_biotype = BiotypeAssignation(self.path_to_output_folder)
		get_biotype.get_information_from_BED_intersection()
		get_biotype.get_biotype_from_intersected_transcripts()
		
		info_peptide_alignments = self.get_info_peptide_alignments()
		self.prepare_info_to_draw_biotypes(info_peptide_alignments, get_biotype.information_final_biotypes_peptides)
	


	def prepare_info_to_draw_biotypes(self, info_peptide_alignments, information_final_biotypes_peptides):

		biotypes_peptides = {}
		data = []
		for type_peptide, peptides in self.input_file_treatment.peptides_by_type.items():
			biotypes_peptides[type_peptide] = {'gene_level_biotype':{}, 'transcript_level_biotype':{}, 'genomic_position_biotype':{}}
			for peptide in peptides:
				try:
					alignments = info_peptide_alignments[peptide]
					count = 0
					for alignment in alignments:
						key = peptide+'_'+alignment[0]
						count = alignment[-1]
						alignment = alignment[0]
						
						try:
							transcripts_intersected = information_final_biotypes_peptides[key]
							for position, transcript in transcripts_intersected.items():
								gene_level_biotype = transcript[0]
								transcript_level_biotype = transcript[1]
								genomic_position_biotype = transcript[2]
								to_add = [type_peptide, peptide, alignment, position, gene_level_biotype, transcript_level_biotype, genomic_position_biotype, count]
						except KeyError:
							biotype_type = 'Intergenic'
							gene_level_biotype = biotype_type
							transcript_level_biotype = biotype_type
							genomic_position_biotype = biotype_type
							to_add = [type_peptide, peptide, alignment, 'No Annotation', gene_level_biotype, transcript_level_biotype, genomic_position_biotype, count]
						
						data.append(to_add)
						if count > 0:
							try:
								biotypes_peptides[type_peptide]['gene_level_biotype'][gene_level_biotype] += 1
							except KeyError:
								biotypes_peptides[type_peptide]['gene_level_biotype'][gene_level_biotype] = 1

							try:
								biotypes_peptides[type_peptide]['transcript_level_biotype'][transcript_level_biotype] += 1
							except KeyError:
								biotypes_peptides[type_peptide]['transcript_level_biotype'][transcript_level_biotype] = 1

							try:
								biotypes_peptides[type_peptide]['genomic_position_biotype'][genomic_position_biotype] += 1
							except KeyError:
								biotypes_peptides[type_peptide]['genomic_position_biotype'][genomic_position_biotype] = 1

				except KeyError:
					pass
		self.draw_biotypes(biotypes_peptides)
		df = pd.DataFrame(data, columns=['Peptide Type', 'Peptide','Alignment','Transcript', 'gene_level_biotype', 'transcript_level_biotype', 'genomic_position_biotype', 'reads count'])
		df.to_csv(self.path_to_output_folder+'/res/annotation_biotypes.csv', index=False, header=True)


	def draw_biotypes(self, biotypes_peptides):

		#others = 'retained_intron', 'IG_V_gene', 'TEC', 'bidirectional_promoter_lncRNA', 'transcribed_unitary_pseudogene', 'transcribed_unprocessed_pseudogene', 'sense_overlapping','processed_pseudogene', 'unprocessed_pseudogene'
		organisation_labels = {'Protein-coding genes':['5UTR', '3UTR', 'In_frame', 'Frameshift', 'protein_coding', 'CDS'], 
							'Non-coding genes':['processed_transcript', 'nonsense_mediated_decay', 'antisense', 'Exons', 'lincRNA', 'Others'], 
							'Intergenic':['Intergenic'], 
							'Intronic':['Introns']}

		for type_peptide, level_biotypes in biotypes_peptides.items():

			for level_biotype, info_level_biotype in level_biotypes.items():

				labels_in_type_peptide = {'Protein-coding genes':{}, 'Non-coding genes': {}, 'Intergenic':{}, 'Intronic':{}}
				outer_labels = []
				outer_sizes = []
				intra_labels = []
				intra_sizes = []
				
				if len(info_level_biotype.values()) > 0:
					title = level_biotype+' '+type_peptide+' peptides'

					for biotype, total_biotype in info_level_biotype.items():
						in_ = False
						for type_, types in organisation_labels.items():
							if biotype in types:
								labels_in_type_peptide[type_][biotype] = total_biotype 
								in_ = True
								break
						if not in_:
							labels_in_type_peptide['Non-coding genes']['Others'] = total_biotype 

					for outer_label_type, intra_labels_type in organisation_labels.items():
						values_in = 0
						for intra_label in intra_labels_type:
							try:
								value = labels_in_type_peptide[outer_label_type][intra_label]
								intra_labels.append(intra_label)
								intra_sizes.append(value)
								values_in += value
							except:
								pass
						if values_in > 0:
							outer_labels.append(outer_label_type)
							outer_sizes.append(values_in)
					
					plots.plot_pie(title, outer_labels, intra_labels, intra_sizes, outer_sizes,  self.path_to_output_folder, self.name_exp, type_peptide+'_'+level_biotype)


	def get_info_peptide_alignments(self):
		info_peptide_alignments = {}
		for peptide_alignment in self.perfect_alignments:
			peptide = peptide_alignment.split('_')[0]
			alignment = peptide_alignment.split('_')[1]
			count = self.perfect_alignments[peptide_alignment][-1]
			try:
				info_peptide_alignments[peptide].append((alignment, count))
			except KeyError:
				info_peptide_alignments[peptide] = [(alignment, count)]
		return info_peptide_alignments




def main(argv):

	path_to_input_folder = ''
	path_to_output_folder = ''
	name_exp = ''
	mode = 'normal'

	try:
		opts, args = getopt.getopt(argv,"hi:n:m:",["path_to_input_folder=", "name_exp=", "mode="])
	except getopt.GetoptError:
		print ('BamQuery.py -i <path_to_input_folder> -n <name_experience> -m <normal/filter>')
		sys.exit(2)

	for opt, arg in opts:
		if opt == '-h':
			print ('BamQuery.py -i <path_to_input_folder> -n <name_experience> -m <normal/filter>')
			sys.exit()
		elif opt in ("-i", "--path_to_input_folder"):
			path_to_input_folder = arg
		elif opt in ("-n", "--name_exp"):
			name_exp = arg
		elif opt in ("-m", "--mode"):
			mode = arg
			if mode == '':
				mode = 'normal'
			elif mode != 'normal' and mode != 'filter':
				print ('BamQuery.py -i <path_to_input_folder> -n <name_experience> -m <normal/filter>')
				sys.exit()
	
	print ('Running BamQuery in mode ', mode)
	if path_to_input_folder[-1] != '/':
		path_to_input_folder += '/'

	path_to_output_folder = '/'.join(path_to_input_folder.split('/')[:-2])+'/output/'
	
	try:
		os.mkdir(path_to_output_folder)
	except OSError:
		print ("The output directory already exist in this path. BamQuery analysis will continue where it left." )
	else:
		print ("Successfully created the directory ", path_to_output_folder)
	
	path_to_logs_folder = path_to_output_folder+'logs/'
	try:
		os.mkdir(path_to_logs_folder)
	except OSError:
		pass

	nameLog = path_to_logs_folder+'BamQuery_Res_'+name_exp+'.log'
	logging.basicConfig(filename=nameLog, filemode='w', level=logging.INFO, format='%(asctime)s %(message)s')
	logging.info('Running BamQuery on experience %s on mode %s ', name_exp, mode)

	path_to_plots_folder = path_to_output_folder+'plots/'
	try:
		os.mkdir(path_to_plots_folder)
	except OSError:
		logging.info("The %s directory already exists. ", path_to_plots_folder)
	else:
		logging.info("Successfully created the directory %s " , path_to_plots_folder)

	path_to_plots_heat_maps_folder = path_to_output_folder+'plots/heat_maps/'
	try:
		os.mkdir(path_to_plots_heat_maps_folder)
	except OSError:
		logging.info("The %s directory already exists. ", path_to_plots_heat_maps_folder)
	else:
		logging.info("Successfully created the directory %s " , path_to_plots_heat_maps_folder)

	path_to_plots_biotypes_folder = path_to_output_folder+'plots/biotypes/'
	try:
		os.mkdir(path_to_plots_biotypes_folder)
	except OSError:
		logging.info("The %s directory already exists. ", path_to_plots_biotypes_folder)
	else:
		logging.info("Successfully created the directory %s " , path_to_plots_biotypes_folder)


	path_to_genome_alignment_folder = path_to_output_folder+'genome_alignments/'
	try:
		os.mkdir(path_to_genome_alignment_folder)
	except OSError:
		logging.info("The %s directory already exists. ", path_to_genome_alignment_folder)
	else:
		logging.info("Successfully created the directory %s " , path_to_genome_alignment_folder)

	path_to_alignment_folder = path_to_output_folder+'alignments/'
	try:
		os.mkdir(path_to_alignment_folder)
	except OSError:
		logging.info("The %s directory already exists. ", path_to_alignment_folder)
	else:
		logging.info("Successfully created the directory %s " , path_to_alignment_folder)

	path_to_res_folder = path_to_output_folder+'res/'
	try:
		os.mkdir(path_to_res_folder)
	except OSError:
		logging.info("The %s directory already exists. ", path_to_res_folder)
	else:
		logging.info("Successfully created the directory %s " , path_to_res_folder)

	path_to_res_bed_files_folder = path_to_output_folder+'res/BED_files/'
	try:
		os.mkdir(path_to_res_bed_files_folder)
	except OSError:
		logging.info("The %s directory already exists. ", path_to_res_bed_files_folder)
	else:
		logging.info("Successfully created the directory %s " , path_to_res_bed_files_folder)


	logging.info('Path to input folder : %s ', path_to_input_folder)
	logging.info('Path to output folder : %s ', path_to_output_folder)
	logging.info('Path to plots folder : %s ', path_to_plots_folder)
	logging.info('Path to genomic alignments STAR folder : %s ', path_to_genome_alignment_folder)
	logging.info('Path to alignments RES folder : %s ', path_to_alignment_folder)
	logging.info('Path to res folder : %s ', path_to_res_folder)

	t0 = time.time()

	BamQuery_obj = BamQuery(path_to_input_folder, path_to_output_folder, name_exp, mode)
	
	t2 = time.time()
	total = t2-t0
	
	logging.info('Total time run function BamQuery to end : %f min', (total/60.0))


if __name__ == "__main__":
	main(sys.argv[1:])


#python /u/ruizma/BAM_Query/Scripts/Python/BamQuery_class.py -p /u/ruizma/BAM_Query/Res/BamQuery_exp_class/peptides.tsv -b /u/ruizma/BAM_Query/Res/BamQuery_exp_class/BAM_directories.tsv -s /u/ruizma/BAM_Query/Res/BamQuery_exp_class/ -n First_Exp
# time python /u/ruizma/BAM_Query/Scripts/Python/BamQuery_class.py -p /u/ruizma/BAM_Query/Res/BamQuery_exp_class/peptides.tsv -b /u/ruizma/BAM_Query/Res/BamQuery_exp_class/BAM_directories.tsv -s /u/ruizma/BAM_Query/Res/BamQuery_exp_class/ -n First_Exp > /u/ruizma/BAM_Query/Res/BamQuery_exp_class/test_reads.reads
# time python /u/ruizma/BAM_Query/Scripts/Python/BamQuery_class.py -p /u/ruizma/BAM_Query/Res/BamQuery_exp_class/peptides.tsv -b /u/ruizma/BAM_Query/Res/BamQuery_exp_class/BAM_directories.tsv -s /u/ruizma/BAM_Query/Res/BamQuery_exp_class/ -n First_Exp
#jupyter3
#time python /u/ruizma/BAM_Query/Scripts/Python/BamQuery/BamQuery.py -i /u/ruizma/BAM_Query/Res/BamQuery_exp_class/Input  -n First_Exp


# time /u/ruizma/BAM_Query/BamQuery/BamQuery_normalisations_5.sh -x full -o /u/ruizma/BAM_Query/Res/BamQuery_exp/ -n test