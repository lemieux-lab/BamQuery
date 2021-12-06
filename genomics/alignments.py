import os,time, subprocess, pickle, multiprocessing, os, csv
import genomics.get_alignments as get_alig
import collections
#from pathos.multiprocessing import ProcessPool
import pandas as pd
import billiard as mp

__author__ = "Maria Virginia Ruiz Cuevas"
__email__ = "maria.virginia.ruiz.cuevas@umontreal.ca"

path_to_lib = '/'.join(os.path.abspath(__file__).split('/')[:-3])+'/lib/'
NUM_WORKERS =  multiprocessing.cpu_count()


def alignment_cs_to_genome(set_peptides, path_to_output_folder, name_exp, light, dbSNP, common, super_logger, var, maxmm, genome_version):

	path_to_output_folder_genome_alignments = path_to_output_folder+'genome_alignments/'
	path_to_output_folder_alignments = path_to_output_folder+'alignments/'

	if genome_version == 'v26_88': 
		index_genome = path_to_lib+'genome_versions/genome_v26_88/Index_BAM_Query/'
	elif genome_version == 'v33_99':
		index_genome = path_to_lib+'genome_versions/genome_v33_99/Index_BAM_Query/'
	else:
		index_genome = path_to_lib+'genome_versions/genome_v38_104/Index_BAM_Query/'

	exist = os.path.exists(path_to_output_folder_genome_alignments+'/Aligned.out.sam')
	exists_light = os.path.exists(path_to_output_folder_alignments+'/Alignments_information_light.dic')
	exists_normal= os.path.exists(path_to_output_folder_alignments+'/Alignments_information.dic')
	exists_fastq= os.path.exists(path_to_output_folder_genome_alignments+'/'+name_exp+'.fastq')

	if not exist and not exists_light and not exists_normal and exists_fastq:
		t_0 = time.time()
		inputFilesR1_1 = path_to_output_folder_genome_alignments+name_exp+'.fastq'
		genomeDirectory = index_genome
		outputFilterMatchInt = 20
		seed = 20

		if maxmm:
			maxMulti = 100000 
			alignTranscriptsPerReadNmax = 30000
			seedPerWindowNmax = 1000
			seedNoneLociPerWindow = 1000
			alignWindowsPerReadNmax = 20000
			alignTranscriptsPerWindowNmax = 1000
		else:
			maxMulti = 10000
			alignTranscriptsPerReadNmax = maxMulti
			seedPerWindowNmax = 500
			seedNoneLociPerWindow = 500
			alignWindowsPerReadNmax = maxMulti
			alignTranscriptsPerWindowNmax = 500

		anchor = maxMulti * 2										
		limitOutSAMoneReadBytes = 33 * 2 * maxMulti
		# https://github.com/alexdobin/STAR/issues/169
		# https://github.com/mhammell-laboratory/TEtranscripts/issues/69
		# https://github.com/alexdobin/STAR/issues/506

		command = 'module add star/2.7.1a; STAR --runThreadN 16'+\
				' --genomeDir '+genomeDirectory+' --seedSearchStartLmax '+str(seed)+\
				' --alignEndsType EndToEnd --sjdbOverhang 32 --sjdbScore 2 --alignSJDBoverhangMin 1 --alignSJoverhangMin 20 '+\
				' --outFilterMismatchNmax 4 --winAnchorMultimapNmax ' +str(anchor)+' --outFilterMultimapNmax '+str(maxMulti)+\
				' --outFilterMatchNmin '+str(outputFilterMatchInt)+ ' --readFilesIn '+inputFilesR1_1+' --outSAMattributes NH HI MD '+\
				' --limitOutSAMoneReadBytes '+str(limitOutSAMoneReadBytes) +\
				' --alignTranscriptsPerWindowNmax '+str(alignTranscriptsPerWindowNmax) +' --alignWindowsPerReadNmax '+str(alignWindowsPerReadNmax)+ \
				' --seedNoneLociPerWindow '+ str(seedNoneLociPerWindow) +' --seedPerWindowNmax '+ str(seedPerWindowNmax)  +\
				' --alignTranscriptsPerReadNmax '+ str(alignTranscriptsPerReadNmax) +' --alignIntronMax 1 --outFileNamePrefix '+path_to_output_folder_genome_alignments
		
		super_logger.info('Command to Align using STAR : %s ', command)
		p_1 = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
		out, err = p_1.communicate()
		t_2 = time.time()
		total = t_2-t_0

		print ("Total time run function alignment_CS_to_genome End : %s min" % (total/60.0))
		super_logger.info('Total time run function alignment_CS_to_genome to end : %f min', (total/60.0))
	else:
		super_logger.info('Alignment file already exists in the output folder : %s --> Skipping this step!', path_to_output_folder_alignments+'/Aligned.out.sam')
	
	perfect_alignments = get_alignments(set_peptides, path_to_output_folder_genome_alignments, path_to_output_folder_alignments, name_exp, light, dbSNP, common, super_logger, var, genome_version)
	
	return perfect_alignments

def get_alignments(set_peptides, path_to_output_folder_genome_alignments, path_to_output_folder_alignments, name_exp, light, dbSNP, common, super_logger, var, genome_version):
	t_0 = time.time()
	sam_file = path_to_output_folder_genome_alignments+'/Aligned.out.sam'
	perfect_alignments = {}
	peptides_with_alignments = set()

	exists_sam_file = os.path.exists(sam_file)
	if exists_sam_file:

		exists = os.path.exists(path_to_output_folder_alignments+'/Alignments_information.dic')
		exists_light = os.path.exists(path_to_output_folder_alignments+'/Alignments_information_light.dic')

		
		if not light and exists_light and not exists:
			perfect_alignments, peptides_with_alignments = filter_peptides_from_alignments_information_light(set_peptides, path_to_output_folder_alignments+'/Alignments_information_light.dic', path_to_output_folder_alignments)
		
		if not exists_light and not exists:

			res_star = get_alig.get_alignments(sam_file, dbSNP, common, super_logger, var, genome_version)
			
			t_2 = time.time()
			total = t_2-t_0

			print ("Total time run function get_alignments End : %s " % (total/60.0))
			super_logger.info('Total time run function get_alignments to end : %f min', (total/60.0))
			super_logger.info('Total perfect aligments : %s ', str(len(res_star[0])))
			
			perfect_alignments = res_star[0]
			variants_alignments = res_star[1]
			
			peptides_with_alignments = res_star[2]
			
			columns = ["Peptide", "Strand", "Alignment", "MCS", "Peptide in Reference", "Diff AA", "Diff ntd", "SNVs"]
			columns_cosmic = ["Peptide", "Strand", "Alignment", "SNV", 'Mutation genome position', 'GRCh', 'Gene', 'SNP', 'Mutation Id', 
							'Mutation CDS', 'Mutation AA', 'Description',
							'Mutation Strand', 'Resistance', 'Score', 'Prediction', 'Status' ]

			alignments = [perfect_alignments, variants_alignments]
			pool = mp.Pool(NUM_WORKERS)
			#pool = ProcessPool(nodes = NUM_WORKERS)
			results = pool.map(generer_alignments_information, alignments)
			
			perfect_alignments_to_print = results[0][0]
			variants_alignments_to_print = results[1][0]
			pool.close()
			pool.join()
			
			info_cosmic = results[0][1]
			
			info_cosmic_to_print = get_info_cosmic(info_cosmic)
			
			if len(info_cosmic_to_print) == 0:
				info_cosmic_to_print = [['NA']*len(columns_cosmic)]

			df1 = pd.DataFrame(perfect_alignments_to_print, columns = columns)
			df2 = pd.DataFrame(variants_alignments_to_print, columns = columns)
			df3 = pd.DataFrame(info_cosmic_to_print, columns = columns_cosmic)

			write_xls_with_alignments_info(path_to_output_folder_alignments, name_exp, df1, df2, df3)
			
			if not light:
				name_path = path_to_output_folder_alignments+'/Alignments_information.dic'
			else:
				name_path = path_to_output_folder_alignments+'/Alignments_information_light.dic'

			with open(name_path, 'wb') as handle:
				pickle.dump(perfect_alignments, handle, protocol=pickle.HIGHEST_PROTOCOL)

			super_logger.info('Alignments Information save to : %s ', name_path)
			
		else:
			super_logger.info('Alignment information already collected in the output folder : %s --> Skipping this step!', path_to_output_folder_alignments+'/Alignments_information.dic')
			
			if not light:
				name_path = path_to_output_folder_alignments+'/Alignments_information.dic'
			else:
				name_path = path_to_output_folder_alignments+'/Alignments_information_light.dic'

			with open(name_path, 'rb') as fp:
				try:
					perfect_alignments = pickle.load(fp)
				except ValueError:
					import pickle5
					perfect_alignments = pickle5.load(fp)

			super_logger.info('Total perfect aligments : %s ', str(len(perfect_alignments)))
			
			try:
				with open(path_to_output_folder_alignments+'alignments/missed_peptides.info') as f:
					for index, line in enumerate(f):
						peptide = line.strip()
						peptides_with_alignments.add(peptide)
			except :
				peptides_keys = list(perfect_alignments.keys())
				for key in peptides_keys:
					peptides_with_alignments.add(key.split('_')[0])

	return perfect_alignments, peptides_with_alignments


def filter_peptides_from_alignments_information_light(peptides, path_to_alignments, path_to_output_folder_alignments):

	with open(path_to_alignments, 'rb') as fp:
		perfect_alignments_light = pickle.load(fp)

	perfect_alignments = {}
	peptides_with_alignments = set()

	for peptide_key, peptide_position_info in perfect_alignments_light.items():
		peptide = peptide_key.split('_')[0]
		if peptide in peptides:
			peptides_with_alignments.add(peptide)
			perfect_alignments[peptide_key] = peptide_position_info

	name_path = path_to_output_folder_alignments+'/Alignments_information.dic'
	with open(name_path, 'wb') as handle:
		pickle.dump(perfect_alignments, handle, protocol=pickle.HIGHEST_PROTOCOL)

	return perfect_alignments, peptides_with_alignments

def generer_alignments_information(alignments_input):

	row_list = []
	info_cosmic = {}

	alignments = collections.OrderedDict(sorted(alignments_input.items()))
	
	for peptide_info, info_alignment in alignments.items():
		alignment = peptide_info.split('_')[1]
		MCS = peptide_info.split('_')[2]
		peptide = peptide_info.split('_')[0]
		
		strand = info_alignment[0]
		peptide_with_snps_local_reference = info_alignment[1]
		
		dif_aas = info_alignment[2]
		snvs = info_alignment[3]
		dif_ntds = info_alignment[4]
		
		snvs_write = ''
		dif_aa_write = ''
		mutations_write = ''
		dif_ntd_write = ''
		
		for snv in snvs:
			snv_to_print = '['+snv[0]+'->'+snv[1]+'|snp:'+snv[2]+'|GenPos:'+str(snv[3])+'|MCSPos:'+str(snv[4])+']'
			snvs_write += snv_to_print
			pos_to_search_cosmic = alignment.split('chr')[1].split(':')[0]+':'+str(snv[3])+'-'+str(snv[3])
			
			try:
				info_cosmic[pos_to_search_cosmic] = [peptide, strand, alignment, snv[2], pos_to_search_cosmic]
			except :
				pass

		for dif_aa in dif_aas:
			dif_aa_write += '['+dif_aa+']'

		for dif_ntd in dif_ntds:
			dif_ntd_write += '['+dif_ntd+']'

		row_list.append([peptide, strand, alignment, MCS, peptide_with_snps_local_reference, dif_aa_write, dif_ntd_write, snvs_write ])

	if len(row_list) == 1048575 :
		return row_list, info_cosmic

	return row_list, info_cosmic


def get_info_cosmic(snv_alignments):

	cosmic_dic = path_to_lib + 'Cosmic_info.dic'

	with open(cosmic_dic, 'rb') as fp:
		cosmic_dic = pickle.load(fp)

	info_cosmic = []

	for positions, peptide_info_snv in snv_alignments.items():
		pos_to_search_cosmic = peptide_info_snv[4]
		try:			
			info_cosmic_ntd = cosmic_dic[pos_to_search_cosmic]
			peptide_info_snv.extend(info_cosmic_ntd)
			info_cosmic_to_print = peptide_info_snv
			info_cosmic.append(info_cosmic_to_print)
		except :
			pass

	return info_cosmic


def write_xls_with_alignments_info(path_to_output_folder_alignments, name_exp, df1, df2, df3):
	writer = pd.ExcelWriter(path_to_output_folder_alignments+name_exp+'_info_alignments.xlsx', engine='xlsxwriter')
	df1.to_excel(writer, sheet_name='Perfect Alignments')
	if len(df2) > 0:
		df2.to_excel(writer, sheet_name='Variants Alignments')
	df3.to_excel(writer, sheet_name='COSMIC Information')
	writer.save()

