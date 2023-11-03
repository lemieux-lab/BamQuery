import os,time, subprocess, pickle, os
import genomics.get_alignments as get_alig
import collections
import pandas as pd
import billiard as mp
from pathos.multiprocessing import ProcessPool
import re, sys, inspect


__author__ = "Maria Virginia Ruiz Cuevas"
__email__ = "maria.virginia.ruiz.cuevas@umontreal.ca"

path_to_lib = '/'.join(os.path.abspath(__file__).split('/')[:-3])+'/lib/'


class NoTraceBackWithLineNumber(Exception):
    def __init__(self, msg):
        try:
            ln = sys.exc_info()[-1].tb_lineno
        except AttributeError:
            ln = inspect.currentframe().f_back.f_lineno
        self.args = "{0.__name__} (line {1}): {2}".format(type(self), ln, msg),
        sys.exit(self)

class NeedMoreInfo(NoTraceBackWithLineNumber):
    pass

def alignment_cs_to_genome(set_peptides, path_to_output_folder, name_exp, light, dbSNP, common, super_logger, var, maxmm, genome_version, mode, mouse, threads):
	global splice_junctions_annotated
	path_to_output_folder_genome_alignments = path_to_output_folder+'genome_alignments/'
	path_to_output_folder_alignments = path_to_output_folder+'alignments/'

	if mouse:
		if genome_version == 'M24':
			index_genome = path_to_lib+'genome_versions/genome_mouse_m24/Index_STAR_2.7.9a/'
			splice_junctions  = index_genome + 'sjdbList.fromGTF.out.tab'

		if genome_version == 'M30':
			index_genome = path_to_lib+'genome_versions/genome_mouse_m30/Index_STAR_2.7.9a/'
			splice_junctions  = index_genome + 'sjdbList.fromGTF.out.tab'
	else:
		if genome_version == 'v26_88': 
			index_genome = path_to_lib+'genome_versions/genome_v26_88/Index_STAR_2.7.9a/'
			splice_junctions  = index_genome + 'sjdbList.fromGTF.out.tab'

		elif genome_version == 'v33_99':
			index_genome = path_to_lib+'genome_versions/genome_v33_99/Index_STAR_2.7.9a/'
			splice_junctions  = index_genome + 'sjdbList.fromGTF.out.tab'

		else:
			index_genome = path_to_lib+'genome_versions/genome_v38_104/Index_STAR_2.7.9a/'
			splice_junctions  = index_genome + 'sjdbList.fromGTF.out.tab'

	splice_junctions_annotated = pd.read_csv(splice_junctions, header=None, sep='\t')
	exist = os.path.exists(path_to_output_folder_genome_alignments+'/Aligned.out.sam')
	exist_sam_dic = os.path.exists(path_to_output_folder_genome_alignments+'/Aligned.out.sam.dic')
	exists_light = os.path.exists(path_to_output_folder_alignments+'/Alignments_information_light.dic')
	exists_normal= os.path.exists(path_to_output_folder_alignments+'/Alignments_information.dic')
	exists_fastq = os.path.exists(path_to_output_folder_genome_alignments+'/'+name_exp+'.fastq')
	
	if not exist and not exist_sam_dic and not exists_light and not exists_normal and exists_fastq:
		t_0 = time.time()
		inputFilesR1_1 = path_to_output_folder_genome_alignments+name_exp+'.fastq'
		genomeDirectory = index_genome
		seed = 20
		maxMulti = 1000000
		anchor = 1000

		if maxmm:
			anchor = 20000
			alignTranscriptsPerReadNmax = 100000
			seedPerWindowNmax = 1500
			seedNoneLociPerWindow = 1500
			alignWindowsPerReadNmax = 20000
			alignTranscriptsPerWindowNmax = 1500
			outFilterMultimapScoreRange = 3
		else:
			alignTranscriptsPerReadNmax = 20000
			seedPerWindowNmax = 1000
			seedNoneLociPerWindow = 1000
			alignWindowsPerReadNmax = 15000
			alignTranscriptsPerWindowNmax = 1000
			outFilterMultimapScoreRange = 2

		if var:
			outFilterMismatchNmax = 4
		else:
			outFilterMismatchNmax = 3

		limitOutSAMoneReadBytes = 2 * ( 33 + 100 ) * maxMulti

		command = 'ulimit -s 8192; STAR --runThreadN '+str(threads)+\
				' --genomeDir '+genomeDirectory+' --seedSearchStartLmax '+str(seed)+\
				' --alignEndsType EndToEnd --sjdbOverhang 32 --alignSJDBoverhangMin 1 --alignSJoverhangMin 10000'+\
				' --outFilterMismatchNmax '+str(outFilterMismatchNmax)+' --outFilterIntronMotifs RemoveNoncanonicalUnannotated --scoreGapNoncan -50 --outFilterType BySJout --winAnchorMultimapNmax '+str(anchor)+' --outFilterMultimapNmax '+str(maxMulti)+\
				' --readFilesIn '+inputFilesR1_1+' --outSAMattributes NH HI MD --limitOutSJcollapsed 5000000'+\
				' --limitOutSAMoneReadBytes '+str(limitOutSAMoneReadBytes) +\
				' --outFilterMultimapScoreRange '+str(outFilterMultimapScoreRange)+' --alignTranscriptsPerWindowNmax '+str(alignTranscriptsPerWindowNmax) +' --alignWindowsPerReadNmax '+str(alignWindowsPerReadNmax)+ \
				' --seedNoneLociPerWindow '+ str(seedNoneLociPerWindow) +' --seedPerWindowNmax '+ str(seedPerWindowNmax)  +\
				' --alignTranscriptsPerReadNmax '+ str(alignTranscriptsPerReadNmax) +' --outFileNamePrefix '+path_to_output_folder_genome_alignments
		
		super_logger.info('Command to Align using STAR : %s ', command)
		p_1 = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
		out, err = p_1.communicate()
		if 'successfully' not in out.decode():
			super_logger.info('A problem occurred while running STAR. \nPlease also verify there is enough memory for the %s number of threads.',str(threads))
			message = '\nA problem occurred while running STAR.\nPlease also verify there is enough memory for the '+ str(threads)+' number of threads.'
			raise NeedMoreInfo(message)
		t_2 = time.time()
		total = t_2-t_0

		print ("Total time run function alignment_CS_to_genome End : %s min" % (total/60.0))
		super_logger.info('Total time run function alignment_CS_to_genome to end : %f min', (total/60.0))
	else:
		super_logger.info('Alignment file already exists in the output folder : %s --> Skipping this step!', path_to_output_folder_genome_alignments+'/Aligned.out.sam')
	
	perfect_alignments = get_alignments(set_peptides, path_to_output_folder_genome_alignments, path_to_output_folder_alignments, name_exp, light, dbSNP, common, super_logger, var, genome_version, mode, mouse, threads)
	return perfect_alignments

def get_alignments(set_peptides, path_to_output_folder_genome_alignments, path_to_output_folder_alignments, name_exp, light, dbSNP, common, super_logger, var, genome_version, mode, mouse, threads):
	t_0 = time.time()
	sam_file = path_to_output_folder_genome_alignments+'/Aligned.out.sam'
	exists = os.path.exists(path_to_output_folder_alignments+'/Alignments_information.dic')
	exists_light = os.path.exists(path_to_output_folder_alignments+'/Alignments_information_light.dic')
		
	perfect_alignments = {}
	peptides_with_alignments = set()

	exists_sam_file = os.path.exists(sam_file)

	if exists_sam_file:

		if not light and exists_light and not exists:
			perfect_alignments, peptides_with_alignments = filter_peptides_from_alignments_information_light(set_peptides, path_to_output_folder_alignments+'/Alignments_information_light_rna.dic', path_to_output_folder_alignments)
		
		if not exists_light and not exists:

			res_star = get_alig.get_alignments(sam_file, path_to_output_folder_genome_alignments, dbSNP, common, super_logger, var, genome_version, mode, mouse, threads)
			
			t_2 = time.time()
			total = t_2-t_0

			print ("Total time run function get_alignments End : %s " % (total/60.0))
			super_logger.info('Total time run function get_alignments to end : %f min', (total/60.0))
			super_logger.info('Total alignments : %s ', str(len(res_star[0])))
			
			perfect_alignments = res_star[0]
			peptides_with_alignments = res_star[1]
			
			columns = ["Peptide", "Strand", "Alignment", "Known Splice Junction", "MCS", "Peptide in Reference", "Diff AA", "Diff ntd", "SNVs"]
			columns_cosmic = ["Peptide", "Strand", "Alignment", "SNV", 'Mutation genome position', 'GRCh', 'Gene', 'SNP', 'Mutation Id', 
							'Mutation CDS', 'Mutation AA', 'Description',
							'Mutation Strand', 'Resistance', 'Score', 'Prediction', 'Status' ]

			pool = ProcessPool(nodes=threads)
			keys = perfect_alignments.keys()
			values = perfect_alignments.values()
			results = pool.map(generation_alignments_information, keys, values)
			perfect_alignments_to_print = []
			info_cosmic = {}
			for res in results:
				perfect_alignments_to_print.extend([res[0]])
				info_cosmic.update(res[1])

			pool.close()
			pool.clear()

			if not mouse and len(info_cosmic) != 0:
				info_cosmic_to_print = get_info_cosmic(info_cosmic)
				
				if len(info_cosmic_to_print) == 0:
					info_cosmic_to_print = [['NA']*len(columns_cosmic)]

				df3 = pd.DataFrame(info_cosmic_to_print, columns = columns_cosmic)
			else:
				df3 = pd.DataFrame()   

			df1 = pd.DataFrame(perfect_alignments_to_print, columns = columns)
			write_xls_with_alignments_info(path_to_output_folder_alignments, name_exp, df1, df3)
			
			if not light:
				name_path = path_to_output_folder_alignments+'Alignments_information.dic'
			else:
				name_path = path_to_output_folder_alignments+'Alignments_information_light.dic'

			with open(name_path, 'wb') as handle:
				pickle.dump(perfect_alignments, handle, protocol=pickle.HIGHEST_PROTOCOL)

			super_logger.info('Alignments Information save to : %s ', name_path)
			
			cols = columns[6:]
			df = df1.drop(cols, axis=1)
			df_to_keep = df.groupby(['Peptide', 'Strand', 'Alignment',  'Known Splice Junction', 'MCS', 'Peptide in Reference']).count().reset_index()

			path = path_to_output_folder_alignments+'/alignments_summary_information.pkl'
			df_to_keep.to_pickle(path)
			
		else:
			super_logger.info('Alignment information already collected in the output folder : %s --> Skipping this step!', path_to_output_folder_alignments+'/Alignments_information.dic')
			
			if not light:
				name_path = path_to_output_folder_alignments+'Alignments_information.dic'
			else:
				name_path = path_to_output_folder_alignments+'Alignments_information_light.dic'


			with open(name_path, 'rb') as fp:
				perfect_alignments = pickle.load(fp)
			
			super_logger.info('Total alignments : %s ', str(len(perfect_alignments)))
			
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


def generation_alignments_information(peptide_info, info_alignment):

	row_info = []
	info_cosmic = {}

	alignment = peptide_info.split('_')[1]
	peptide = peptide_info.split('_')[0]

	try:
		MCS = peptide_info.split('_')[2]
	except IndexError:
		MCS = ''
	strand = info_alignment[0]
	peptide_with_snps_local_reference = info_alignment[1]
	
	dif_aas = info_alignment[2]
	snvs = info_alignment[3]
	dif_ntds = info_alignment[4]

	chr = alignment.split(':')[0]
	known_splice_junction = []

	if '|' in alignment:
		result = re.findall(r"\d+", alignment.split(':')[1])[1:-1]
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
	snvs_write = ''
	dif_aa_write = ''
	mutations_write = ''
	dif_ntd_write = ''
	
	for snv in snvs:
		snv_to_print = snv[0]+'>'+snv[1]+'|snp:'+snv[2]+'|GenPos:'+str(snv[3])+'|MCSPos:'+str(snv[4])+','
		snvs_write += snv_to_print
		pos_to_search_cosmic = alignment.split('chr')[1].split(':')[0]+':'+str(snv[3])+'-'+str(snv[3])
		
		try:
			info_cosmic[pos_to_search_cosmic] = [peptide, strand, alignment, snv[2], pos_to_search_cosmic]
		except :
			pass

	for dif_aa in dif_aas:
		dif_aa_write += dif_aa+','

	for dif_ntd in dif_ntds:
		dif_ntd_write += dif_ntd+','

	row_info = [peptide, strand, alignment, known_splice_junction, MCS, peptide_with_snps_local_reference, dif_aa_write[:-1], dif_ntd_write[:-1], snvs_write[:-1] ]
	
	return row_info, info_cosmic


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


def write_xls_with_alignments_info(path_to_output_folder_alignments, name_exp, df1, df3):
	if len(df1) < 1048576:
		writer = pd.ExcelWriter(path_to_output_folder_alignments+name_exp+'_info_alignments.xlsx', engine='xlsxwriter')
		writer.book.use_zip64()
		df1.to_excel(writer, sheet_name='Perfect Alignments', index=None)
		df3.to_excel(writer, sheet_name='COSMIC Information', index=None)
		writer.save()
	else:
		df1.to_csv(path_to_output_folder_alignments+name_exp+'_info_alignments.csv', index=None)
		df3.to_csv(path_to_output_folder_alignments+name_exp+'_cosmic_information.csv', index=None)
		

