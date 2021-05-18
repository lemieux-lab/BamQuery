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
	logging.basicConfig(filename=nameLog, filemode='a', level=logging.INFO, format='%(asctime)s %(message)s')
	logging.info('=============== BamQuery id : %s, mode : %s, strandedness : %s ===================', name_exp, mode, strandedness)

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


	path_to_plots_biotypes_folder = path_to_output_folder+'plots/biotypes/'
	try:
		os.mkdir(path_to_plots_biotypes_folder)
	except OSError:
		pass

	path_to_plots_heat_maps_genome_annotation_folder = path_to_output_folder+'plots/biotypes/genome_annotation/'
	try:
		os.mkdir(path_to_plots_heat_maps_genome_annotation_folder)
	except OSError:
		pass

	path_to_plots_heat_maps_ERE_annotation_folder = path_to_output_folder+'plots/biotypes/ERE_annotation/'
	try:
		os.mkdir(path_to_plots_heat_maps_ERE_annotation_folder)
	except OSError:
		pass

	path_to_plots_heat_maps_genome_ERE_annotation_folder = path_to_output_folder+'plots/biotypes/genome_and_ERE_annotation/'
	try:
		os.mkdir(path_to_plots_heat_maps_genome_ERE_annotation_folder)
	except OSError:
		pass

	path_to_plots_heat_maps_biotype_by_peptide_type_reads_folder = path_to_output_folder+'plots/biotypes/genome_annotation/biotype_by_peptide_type_reads/'
	try:
		os.mkdir(path_to_plots_heat_maps_biotype_by_peptide_type_reads_folder)
	except OSError:
		pass

	path_to_plots_heat_maps_by_level_biotype_folder = path_to_output_folder+'plots/biotypes/genome_annotation/biotype_by_peptide_type_reads/by_level_biotype/'
	try:
		os.mkdir(path_to_plots_heat_maps_by_level_biotype_folder)
	except OSError:
		pass

	path_to_plots_heat_maps_gene_level_biotype_folder = path_to_output_folder+'plots/biotypes/genome_annotation/biotype_by_peptide_type_reads/by_level_biotype/gene_level_biotype/'
	try:
		os.mkdir(path_to_plots_heat_maps_gene_level_biotype_folder)
	except OSError:
		pass

	path_to_plots_heat_maps_transcript_level_biotype_folder = path_to_output_folder+'plots/biotypes/genome_annotation/biotype_by_peptide_type_reads/by_level_biotype/transcript_level_biotype/'
	try:
		os.mkdir(path_to_plots_heat_maps_transcript_level_biotype_folder)
	except OSError:
		pass

	path_to_plots_heat_maps_genomic_position_biotype_folder = path_to_output_folder+'plots/biotypes/genome_annotation/biotype_by_peptide_type_reads/by_level_biotype/genomic_position_biotype/'
	try:
		os.mkdir(path_to_plots_heat_maps_genomic_position_biotype_folder)
	except OSError:
		pass

	path_to_plots_heat_maps_all_peptides_biotype_folder = path_to_output_folder+'plots/biotypes/genome_annotation/biotype_by_peptide_type_reads/all_peptides/'
	try:
		os.mkdir(path_to_plots_heat_maps_all_peptides_biotype_folder)
	except OSError:
		pass

	path_to_plots_heat_maps_global_biotypes_folder = path_to_output_folder+'plots/biotypes/genome_annotation/global_biotypes/'
	try:
		os.mkdir(path_to_plots_heat_maps_global_biotypes_folder)
	except OSError:
		pass

	path_to_plots_heat_maps_global_biotypes_cat_folder = path_to_output_folder+'plots/biotypes/genome_annotation/global_biotypes/by_peptide_type/'
	try:
		os.mkdir(path_to_plots_heat_maps_global_biotypes_cat_folder)
	except OSError:
		pass

	path_to_plots_heat_maps_global_biotypes_all_folder = path_to_output_folder+'plots/biotypes/genome_annotation/global_biotypes/all_peptides/'
	try:
		os.mkdir(path_to_plots_heat_maps_global_biotypes_all_folder)
	except OSError:
		pass

	path_to_plots_heat_maps_genome_ERE_annotation_all_folder = path_to_output_folder+'plots/biotypes/genome_and_ERE_annotation/all_peptides/'
	try:
		os.mkdir(path_to_plots_heat_maps_genome_ERE_annotation_all_folder)
	except OSError:
		pass

	path_to_plots_heat_maps_genome_ERE_annotation_by_peptide_type_folder = path_to_output_folder+'plots/biotypes/genome_and_ERE_annotation/by_peptide_type/'
	try:
		os.mkdir(path_to_plots_heat_maps_genome_ERE_annotation_by_peptide_type_folder)
	except OSError:
		pass

	path_to_plots_heat_maps_ERE_global_biotype_folder = path_to_output_folder+'plots/biotypes/ERE_annotation/global_biotypes/'
	try:
		os.mkdir(path_to_plots_heat_maps_ERE_global_biotype_folder)
	except OSError:
		pass

	path_to_plots_heat_maps_ERE_global_biotypes_biotype_folder = path_to_output_folder+'plots/biotypes/ERE_annotation/with_counts/'
	try:
		os.mkdir(path_to_plots_heat_maps_ERE_global_biotypes_biotype_folder)
	except OSError:
		pass

	path_to_plots_correlation_folder = path_to_output_folder+'plots/correlation/'
	try:
		os.mkdir(path_to_plots_correlation_folder)
	except OSError:
		pass

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

	path_to_res_folder = path_to_output_folder+'res/'
	try:
		os.mkdir(path_to_res_folder)
	except OSError:
		pass

	path_to_res_bed_files_folder = path_to_output_folder+'res/BED_files/'
	try:
		os.mkdir(path_to_res_bed_files_folder)
	except OSError:
		pass

	path_to_res_aux_files_folder = path_to_output_folder+'res/AUX_files/'
	try:
		os.mkdir(path_to_res_aux_files_folder)
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


	logging.info('Path to input folder : %s ', path_to_input_folder)
	logging.info('Path to output folder : %s ', path_to_output_folder)
	logging.info('Path to plots folder : %s ', path_to_plots_folder)
	logging.info('Path to genomic alignments STAR folder : %s ', path_to_genome_alignment_folder)
	logging.info('Path to alignments RES folder : %s ', path_to_alignment_folder)
	logging.info('Path to res folder : %s ', path_to_res_folder)
	logging.info('=============== # ===================')

	return path_to_output_folder
