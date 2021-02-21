import warnings, time, sys, os, datetime, getopt, logging, _thread
warnings.filterwarnings("ignore")
from os import listdir
from os.path import isfile, join

from readers.read_input import ReadInputFile
from readers.info_transcripts_annotation import InfoTranscripts
from readers.intersection_alignments_annotations import IntersectAnnotations

from utils.primary_read_count import GetPrimaryReadCountBamFiles
from utils.reverse_translation import ReverseTranslation

from genomics.alignments import Alignments
from genomics.get_counts import GetCounts
from genomics.normalization import Normalization
from genomics.get_biotype import BiotypeAssignation

import plotting.heat_maps as heat_maps

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
		
		self.perfect_alignments = self.res_star[0]
		get_counts = GetCounts(self.path_to_output_folder, self.name_exp, self.input_file_treatment.manual_mode, self.mode)
		df_counts = get_counts.get_counts(self.perfect_alignments, self.bam_files_info.bam_files_list)
		heat_maps.get_heat_map(df_counts, self.path_to_output_folder, self.name_exp, 'get_counts')

		normalization = Normalization(self.path_to_output_folder, self.name_exp, self.bam_files_info.bam_files_list, self.mode)
		def_norm = normalization.get_normalization(df_counts, '_norm.csv')
		heat_maps.get_heat_map(def_norm, self.path_to_output_folder, self.name_exp, 'get_norm')


	def run_bam_query_filter_mode(self):
		self.common_to_modes()

		get_counts = GetCounts(self.path_to_output_folder, self.name_exp, self.input_file_treatment.manual_mode, self.mode)
		perfect_alignments_to_return, df_counts = get_counts.filter_counts(self.res_star[0], self.bam_files_info.bam_ribo_files_list)
		heat_maps.get_heat_map(df_counts, self.path_to_output_folder, self.name_exp, 'filter_counts')
		
		normalization = Normalization(self.path_to_output_folder, self.name_exp, self.bam_files_info.bam_ribo_files_list, self.mode)
		def_norm = normalization.get_normalization(df_counts, '_filter_norm.csv')
		heat_maps.get_heat_map(def_norm, self.path_to_output_folder, self.name_exp, 'filter_norm')

		self.perfect_alignments = perfect_alignments_to_return

		df_counts = get_counts.get_counts(self.perfect_alignments, self.bam_files_info.bam_files_list)
		heat_maps.get_heat_map(df_counts, self.path_to_output_folder, self.name_exp, '_counts')

		normalization = Normalization(self.path_to_output_folder, self.name_exp, self.bam_files_info.bam_files_list, self.mode)
		def_norm = normalization.get_normalization(df_counts, '_norm.csv')
		heat_maps.get_heat_map(def_norm, self.path_to_output_folder, self.name_exp, '_norm')


	def common_to_modes(self):
		self.bam_files_info = GetPrimaryReadCountBamFiles()
		self.bam_files_info.get_all_counts(self.path_to_input_folder, self.path_to_output_folder)

		self.input_file_treatment = ReadInputFile(self.path_to_input_folder)
		self.input_file_treatment.treatment_file()
		
		self.reverse_translation = ReverseTranslation()
		self.reverse_translation.reverse_translation(self.input_file_treatment.peptide_mode, self.input_file_treatment.CS_mode, self.path_to_output_folder, self.name_exp)
		
		self.alignments = Alignments(self.path_to_output_folder, self.name_exp)
		self.res_star = self.alignments.alignment_cs_to_genome()

	
	def get_annotations(self):
		get_info_transcripts = InfoTranscripts()
		get_info_transcripts.set_values()

		intsersect_to_annotations = IntersectAnnotations(self.perfect_alignments, self.path_to_output_folder, self.name_exp)
		intsersect_to_annotations.generate_BED_files()
		intsersect_to_annotations.perform_intersection_with_annotation()

		get_biotype = BiotypeAssignation(self.path_to_output_folder)
		get_biotype.get_information_from_BED_intersection()

def main(argv):

	path_to_input_folder = ''
	path_to_output_folder = ''
	name_exp = ''
	mode = ''

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
			if mode != 'normal' and mode != 'filter':
				print ('BamQuery.py -i <path_to_input_folder> -n <name_experience> -m <normal/filter>')
				sys.exit()
	
	if path_to_input_folder[-1] != '/':
		path_to_input_folder += '/'

	path_to_output_folder = '/'.join(path_to_input_folder.split('/')[:-2])+'/output/'
	
	try:
		os.mkdir(path_to_output_folder)
	except OSError:
		print ("The %s directory already exist in this path. " % path_to_output_folder)
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