import os, logging

__author__ = "Maria Virginia Ruiz"
__email__ = "maria.virginia.ruiz.cuevas@umontreal.ca"


def directories_creation(path_to_input_folder, name_exp, mode, light):

	path_to_output_folder = '/'.join(path_to_input_folder.split('/')[:-2])+'/output/'
	
	try:
		os.mkdir(path_to_output_folder)
	except OSError:
		print ("The output directory already exists in this path. \nBamQuery analysis will continue where it left." )
	else:
		print ("Successfully created the directory ", path_to_output_folder)
	
	path_to_logs_folder = path_to_output_folder+'logs/'
	try:
		os.mkdir(path_to_logs_folder)
	except OSError:
		pass

	formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')


	def setup_logger(name, log_file, mode, level=logging.INFO):
		
		handler = logging.FileHandler(log_file, mode=mode)        
		handler.setFormatter(formatter)

		logger = logging.getLogger(name)
		logger.setLevel(level)
		logger.addHandler(handler)

		return logger
	
	# First logger
	nameLog = path_to_logs_folder+'BamQuery_Res_'+name_exp+'.log'
	exists = os.path.exists(nameLog)
	
	super_logger = setup_logger('Super_Logger', nameLog, 'a')
	if exists:
		super_logger.info('')
		super_logger.info('')
		super_logger.info('BamQuery analysis will continue where it left....')
	
	# Second logger
	nameLog = path_to_logs_folder+'Information_BAM_directories.log'
	bam_files_logger = setup_logger('Bam_Files_Logger', nameLog, 'a')
	
	if mode == 'translation':
		paths = [path_to_output_folder+'res_translation/temps_files/', 
		path_to_output_folder+'plots/heat_maps/translation_evidence_heatmap/total_translation_expression_heatmap/',
		path_to_output_folder+'plots/heat_maps/translation_evidence_heatmap/average_translation_expression_heatmap/', path_to_output_folder+'res_translation/BED_files/', path_to_output_folder+'res_translation/AUX_files/']
		
		for path in paths:
			try:
				os.makedirs(path)
			except OSError:
				pass

	elif not light :

		paths = [path_to_output_folder+'res/BED_files/', 
		path_to_output_folder+'plots/heat_maps/transcription_evidence_heatmap/total_transcription_expression_heatmap/', 
		path_to_output_folder+'plots/heat_maps/transcription_evidence_heatmap/average_transcription_expression_heatmap/', 
		path_to_output_folder+'plots/biotypes/biotype_by_sample_group/all_peptides/', 
		path_to_output_folder+'plots/biotypes/biotype_by_sample_group/by_peptide_type/',
		path_to_output_folder+'plots/biotypes/genome_and_ERE_annotation/all_peptides/',
		path_to_output_folder+'plots/biotypes/genome_and_ERE_annotation/by_peptide_type/',
		path_to_output_folder+'res/biotype_classification/full_info_biotypes/',
		path_to_output_folder+'res/biotype_classification/summary_info_biotypes/',
		path_to_output_folder+'res/temps_files/']

		for path in paths:
			try:
				os.makedirs(path)
			except OSError:
				pass

	else:
		paths = [path_to_output_folder+'res_light/temps_files/']
		
		for path in paths:
			try:
				os.makedirs(path)
			except OSError:
				pass

	if not os.path.exists(path_to_output_folder+'alignments/Alignments_information_light.dic')	:
		path_to_genome_alignment_folder = path_to_output_folder+'genome_alignments/'
		try:
			os.mkdir(path_to_genome_alignment_folder)
		except OSError:
			pass

	path_to_alignment_folder = path_to_output_folder+'alignments/'
	try:
		os.mkdir(path_to_alignment_folder)
	except OSError:
		pass

	
	
	super_logger.info('Path to input folder : %s ', path_to_input_folder)
	super_logger.info('Path to output folder : %s ', path_to_output_folder)
	super_logger.info('=============== # ===================')

	return path_to_output_folder, super_logger, bam_files_logger


