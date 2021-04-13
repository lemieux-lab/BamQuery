import os, logging

__author__ = "Maria Virginia Ruiz"
__email__ = "maria.virginia.ruiz.cuevas@umontreal.ca"


def directories_creation(path_to_input_folder, name_exp, mode, strandedness):

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

	nameLog = path_to_logs_folder+'BamQuery_Res_'+name_exp+'.log'
	logging.basicConfig(filename=nameLog, filemode='w', level=logging.INFO, format='%(asctime)s %(message)s')
	logging.info('=============== BamQuery id : %s, mode : %s, strandedness : %s ===================', name_exp, mode, strandedness)

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

	path_to_plots_correlation_folder = path_to_output_folder+'plots/correlation/'
	try:
		os.mkdir(path_to_plots_correlation_folder)
	except OSError:
		logging.info("The %s directory already exists. ", path_to_plots_correlation_folder)
	else:
		logging.info("Successfully created the directory %s " , path_to_plots_correlation_folder)


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
	logging.info('=============== # ===================')

	return path_to_output_folder
