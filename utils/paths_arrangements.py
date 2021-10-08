import os, logging

__author__ = "Maria Virginia Ruiz"
__email__ = "maria.virginia.ruiz.cuevas@umontreal.ca"


def directories_creation(path_to_input_folder, name_exp, mode, strandedness, light):

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
	
	
	if not light :

		path_to_res_folder = path_to_output_folder+'res/'
		try:
			os.mkdir(path_to_res_folder)
		except OSError:
			pass

		path_to_plots_folder = path_to_output_folder+'plots/'
		try:
			os.mkdir(path_to_plots_folder)
		except OSError:
			pass
	
		path_to_plots_heat_maps_folder = path_to_output_folder+'plots/heat_maps/'
		try:
			os.mkdir(path_to_plots_heat_maps_folder)
		except OSError:
			pass

		path_to_res_bed_files_folder = path_to_output_folder+'res/BED_files/'
		try:
			os.mkdir(path_to_res_bed_files_folder)
		except OSError:
			pass

		path_to_plots_biotypes_folder = path_to_output_folder+'plots/biotypes/'
		try:
			os.mkdir(path_to_plots_biotypes_folder)
		except OSError:
			pass

		path_to_plots_biotype_by_sample_group_folder = path_to_output_folder+'plots/biotypes/biotype_by_sample_group/'
		try:
			os.mkdir(path_to_plots_biotype_by_sample_group_folder)
		except OSError:
			pass

		path_to_plots_biotype_by_sample_group_all_folder = path_to_output_folder+'plots/biotypes/biotype_by_sample_group/all_peptides/'
		try:
			os.mkdir(path_to_plots_biotype_by_sample_group_all_folder)
		except OSError:
			pass

		path_to_plots_biotype_by_sample_group_by_peptide_type_folder = path_to_output_folder+'plots/biotypes/biotype_by_sample_group/by_peptide_type/'
		try:
			os.mkdir(path_to_plots_biotype_by_sample_group_by_peptide_type_folder)
		except OSError:
			pass
			
		path_to_plots_biotype_genome_ERE_annotation_folder = path_to_output_folder+'plots/biotypes/genome_and_ERE_annotation/'
		try:
			os.mkdir(path_to_plots_biotype_genome_ERE_annotation_folder)
		except OSError:
			pass

		path_to_plots_biotype_genome_ERE_annotation_all_folder = path_to_output_folder+'plots/biotypes/genome_and_ERE_annotation/all_peptides/'
		try:
			os.mkdir(path_to_plots_biotype_genome_ERE_annotation_all_folder)
		except OSError:
			pass

		path_to_plots_biotype_genome_ERE_annotation_by_peptide_type_folder = path_to_output_folder+'plots/biotypes/genome_and_ERE_annotation/by_peptide_type/'
		try:
			os.mkdir(path_to_plots_biotype_genome_ERE_annotation_by_peptide_type_folder)
		except OSError:
			pass

		path_to_res_aux_files_folder = path_to_output_folder+'res/AUX_files/'
		try:
			os.mkdir(path_to_res_aux_files_folder)
		except OSError:
			pass

		path_to_res_full_info_folder = path_to_output_folder+'res/full_info_biotypes/'
		try:
			os.mkdir(path_to_res_full_info_folder)
		except OSError:
			pass

		path_to_res_summary_info_folder = path_to_output_folder+'res/summary_info_biotypes/'
		try:
			os.mkdir(path_to_res_summary_info_folder)
		except OSError:
			pass

		path_to_res_temps_files_folder = path_to_output_folder+'res/temps_files/'
		try:
			os.mkdir(path_to_res_temps_files_folder)
		except OSError:
			pass

		path_to_res_aux_processed_files_folder = path_to_output_folder+'res/AUX_files/processed/'
		try:
			os.mkdir(path_to_res_aux_processed_files_folder)
		except OSError:
			pass

		path_to_res_aux_processed_rna_files_folder = path_to_output_folder+'res/AUX_files/processed/rna_norm/'
		try:
			os.mkdir(path_to_res_aux_processed_rna_files_folder)
		except OSError:
			pass

	else:
		path_to_res_folder = path_to_output_folder+'res_light/'
		try:
			os.mkdir(path_to_res_folder)
		except OSError:
			pass

		path_to_res_bed_files_folder = path_to_output_folder+'res_light/BED_files/'
		try:
			os.mkdir(path_to_res_bed_files_folder)
		except OSError:
			pass
			
		path_to_res_temps_files_folder = path_to_output_folder+'res_light/temps_files/'
		try:
			os.mkdir(path_to_res_temps_files_folder)
		except OSError:
			pass

		path_to_res_aux_files_folder = path_to_output_folder+'res_light/AUX_files/'
		try:
			os.mkdir(path_to_res_aux_files_folder)
		except OSError:
			pass

		path_to_res_aux_processed_files_folder = path_to_output_folder+'res_light/AUX_files/processed/'
		try:
			os.mkdir(path_to_res_aux_processed_files_folder)
		except OSError:
			pass

		path_to_res_aux_processed_rna_files_folder = path_to_output_folder+'res_light/AUX_files/processed/rna_norm/'
		try:
			os.mkdir(path_to_res_aux_processed_rna_files_folder)
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

	
	if mode == 'translation':
		path_to_res_aux_processed_ribo_files_folder = path_to_output_folder+'res/AUX_files/processed/ribo_norm/'
		try:
			os.mkdir(path_to_res_aux_processed_ribo_files_folder)
		except OSError:
			pass

		path_to_res_aux_processed_rna_ribo_files_folder = path_to_output_folder+'res/AUX_files/processed/rna_ribo_norm/'
		try:
			os.mkdir(path_to_res_aux_processed_rna_ribo_files_folder)
		except OSError:
			pass

		path_to_plots_correlation_folder = path_to_output_folder+'plots/correlation/'
		try:
			os.mkdir(path_to_plots_correlation_folder)
		except OSError:
			pass


	super_logger.info('Path to input folder : %s ', path_to_input_folder)
	super_logger.info('Path to output folder : %s ', path_to_output_folder)
	super_logger.info('=============== # ===================')

	return path_to_output_folder, super_logger, bam_files_logger
