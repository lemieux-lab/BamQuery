import warnings
warnings.filterwarnings("ignore")
import time, os
import pysam, re
from itertools import groupby
from operator import itemgetter
import utils.useful_functions as uf
import pickle
from pathos.multiprocessing import ProcessPool
import multiprocessing
import collections

__author__ = "Maria Virginia Ruiz Cuevas"
__email__ = "maria.virginia.ruiz.cuevas@umontreal.ca"

path_to_lib = '/'.join(os.path.abspath(__file__).split('/')[:-3])+'/lib/'

def set_strand_read(strand):
	number = "{0:b}".format(int(strand))
	if len(number)>= 5:
		if '1' == number[-5]:
			return '-'
		else:
			return '+'
	else:
		return '+'

def get_ranges(cigar, start, lenSeq, strand, chr, faFile):

	rang = [0]*lenSeq
	rx = re.findall('(\d+)([MISDNX=])?', cigar)
	realReadStart = start
	realEnd = start
	operators = []
	lastIndex = 0
	seqReference = ''

	for index in rx:
		operation = index[1]
		length = int(index[0])
		operators.append(operation)

		if ('S' in operation):
			end = length
			lastIndex += end
			realReadStart = start - length
			realEnd += end
			resSeq = faFile.fetch(chr, start-1,start+end-1)
			seqReference += resSeq

		elif ('I' in index) or ('M' in index) or ('=' in index) or ('X' in index) :
			end = length
			resSeq = faFile.fetch(chr, start-1, start+end-1)
			seqReference += resSeq
			
			for i in range(0,end):
				rang[lastIndex+i] = start
				start = start + 1

			realEnd += end
			lastIndex = lastIndex + end

		elif ('N' in index) or ('D' in index):
			end = start + length
			start = end
			realEnd += length

	if strand == '-':
		seqReference = uf.reverseComplement(seqReference)
		rang = rang[::-1]

	if 'N' in seqReference:
		seqReference = seqReference.replace("N", "C")
	if 'n' in seqReference:
		seqReference = seqReference.replace("n", "C")
	
	return rang, operators, seqReference


def get_local_reference(start, lenSeq, chr, strand, faFile):

	seqReference = ''
	end = start+lenSeq-1
	rang = range(start, end+1)
	seqReference = faFile.fetch(chr, start-1,end)

	if strand == '-':
		seqReference = uf.reverseComplement(seqReference)
		rang = rang[::-1]

	if 'N' in seqReference:
		seqReference = seqReference.replace("N", "C")
	if 'n' in seqReference:
		seqReference = seqReference.replace("n", "C")
	
	return rang, seqReference


def find_ranges(iterable):
	ranges = []
	ordered = sorted(iterable)
	toReturn = ''
	for k, g in groupby(enumerate(ordered), key=lambda t :t[1] - t[0]):
		group = list(map(itemgetter(1), g))
		ranges.append((group[0], group[-1]))

	for pos in ranges:
		if pos[0] != 0 and pos[1] != 0:
			toReturn += str(pos[0])+'-'+str(pos[1])+'|'

	toReturn = toReturn[:-1]
	return toReturn

def read_sam_file(sam_file):
	time0 = time.time()
	alignments_by_chromosome_strand = {}
	MCS_alignments_by_chromosome_strand = {}
	
	cont = 0
	with open(sam_file) as f:
		for index, line in enumerate(f):
			if '@' not in line:
				splitLine = line.strip().split('\t')
				queryname = splitLine[0]
				cigar = splitLine[5]
				chr = splitLine[2]
				readStart = int(splitLine[3])
				strand = set_strand_read(splitLine[1])
				peptide = queryname.split('_')[0]
				#number_MCS = queryname.split('_')[1]
				number_hit = int(splitLine[11].split('NH')[1].split(':')[2].split('\t')[0])
				MCS = splitLine[9]
				lenSeq = len(MCS)
				
				if chr.isdigit() or chr == 'X' or chr == 'Y' or chr == 'M' or chr == 'MT':
					if chr == 'MT':
						chr = 'chrM'
					else:
						chr = 'chr'+chr

				if strand == '-':
					MCS = uf.reverseComplement(MCS)
				
				key = str(readStart)+'|'+cigar+'|'+MCS+'|'+peptide+'|'+strand#+'|'+number_MCS
				
				try:
					chromosomes_alignments = alignments_by_chromosome_strand[chr]
					chromosomes_alignments.add(key)
				except KeyError:
					set_alignments = set()
					set_alignments.add(key)
					alignments_by_chromosome_strand[chr] = set_alignments
					
	timeFinal = time.time()
	total = (timeFinal-time0) / 60.0
	total_chromosomes = len(alignments_by_chromosome_strand)
	super_logger.info('Time reading SAM file : %f min Chromosomes %d ', total, total_chromosomes)
	name_path = sam_file+'.dic'
	with open(name_path, 'wb') as handle:
		pickle.dump(alignments_by_chromosome_strand, handle, protocol=pickle.HIGHEST_PROTOCOL)

	return alignments_by_chromosome_strand


def get_alignments_chromosome(chr, chromosomes_alignments):

	positions_mcs_peptides_perfect_alignment = {}
	positions_mcs_peptides_variants_alignment = {}
	
	chromosome = {}
	peptides_in = set()

	local_visited = set()
	faFile = pysam.FastaFile(genomePath, genomePathFai)

	if path_to_db != '':
		try:
			with open(path_to_db+chr+'.dic', 'rb') as fp:
				chromosome = pickle.load(fp)
		except IOError:
			chromosome = {}
	else:
		chromosome = {}
		
	for position in chromosomes_alignments: 
		split_position = position.split('|')
		
		readStart = int(split_position[0])
		cigar = split_position[1]
		MCS = split_position[2]
		peptide = split_position[3]
		strand  = split_position[4]
		lenSeq = len(MCS)

		end = readStart+lenSeq-1
		positions_trans = find_ranges(range(readStart, end+1))
		range_trans_local = chr+':'+positions_trans
		key_local = peptide+'_'+range_trans_local
		
		rang, operators, seq_reference_align = get_ranges(cigar, readStart, lenSeq, strand, chr, faFile)
		positions_trans = find_ranges(rang)
		range_trans = chr+':'+positions_trans
		key = peptide+'_'+range_trans
		
		if key_local not in local_visited:
			local_visited.add(key_local)
			rang_local_ref, seq_reference_local = get_local_reference(readStart, lenSeq, chr, strand, faFile)
			MCS_perfect_alignments_local, MCS_variant_alignments_local = get_sequences_at_position_local(peptide, seq_reference_local, MCS, list(rang_local_ref), strand, chromosome, chr)
			
			if len(MCS_perfect_alignments_local) > 0:
				peptides_in.add(peptide)
			
				MCS_in =  MCS_perfect_alignments_local[0]
				key_local = key_local+'_'+MCS_in
				
				local_translation_peptide = MCS_perfect_alignments_local[1][0]
				differences_pep = MCS_perfect_alignments_local[1][1]
				info_snps = MCS_perfect_alignments_local[1][2]
				differences_ntds = MCS_perfect_alignments_local[1][3]

				positions_mcs_peptides_perfect_alignment[key_local] = [strand, local_translation_peptide, differences_pep, info_snps, differences_ntds, []]
				
			if len(MCS_variant_alignments_local) > 0:
				MCS_in_var =  MCS_variant_alignments_local[0]
				key_local = key_local+'_'+MCS_in_var

				local_translation_peptide = MCS_variant_alignments_local[1][0]
				differences_pep = MCS_variant_alignments_local[1][1]
				info_snps = MCS_variant_alignments_local[1][2]
				differences_ntds = MCS_variant_alignments_local[1][3]

				positions_mcs_peptides_variants_alignment[key_local] = [strand, local_translation_peptide, differences_pep, info_snps, differences_ntds, []]

		 
		MCS_perfect_alignments_align, MCS_variant_alignments_align = get_sequences_at_position(peptide, seq_reference_align, MCS, rang, strand, chromosome, chr)
		
		if len(MCS_perfect_alignments_align) > 0:
			peptides_in.add(peptide)

			MCS_in =  MCS_perfect_alignments_align[0]
			key = key+'_'+MCS_in

			local_translation_peptide = MCS_perfect_alignments_align[1][0]
			differences_pep = MCS_perfect_alignments_align[1][1]
			info_snps = MCS_perfect_alignments_align[1][2]
			differences_ntds = MCS_perfect_alignments_align[1][3]

			positions_mcs_peptides_perfect_alignment[key] = [strand, local_translation_peptide, differences_pep, info_snps, differences_ntds, []]
		

		if len(MCS_variant_alignments_align) > 0:
			MCS_in_var =  MCS_variant_alignments_align[0]
			key = key+'_'+MCS

			local_translation_peptide = MCS_variant_alignments_align[1][0]
			differences_pep = MCS_variant_alignments_align[1][1]
			info_snps = MCS_variant_alignments_align[1][2]
			differences_ntds = MCS_variant_alignments_align[1][3]

			positions_mcs_peptides_variants_alignment[key] = [strand, local_translation_peptide, differences_pep, info_snps, differences_ntds, []]

	chromosome = {}
	return positions_mcs_peptides_perfect_alignment, positions_mcs_peptides_variants_alignment, peptides_in


def reverse_translation(ntd):
	if ntd == 'A':
		return 'T'
	elif ntd == 'C':
		return 'G'
	elif ntd == 'G':
		return 'C'
	else:
		return 'A'


def get_sequences_at_position(peptide, seq_reference_local, MCS, rang_local_ref, strand, chromosome, chr):

	MCS_perfect_alignments = []
	MCS_variant_alignments = []
	
	local_translation_peptide = translation_seq(chr, seq_reference_local)
	differences_ntds = [seq_reference_local[i]+':'+str(i) for i in range(len(seq_reference_local)) if seq_reference_local[i]!= MCS[i]]
	list_seq_reference_local = list(seq_reference_local)

	info_snps = []
	snps_set = set()
	
	for dif in differences_ntds:
		dif = int(dif.split(':')[1])
		pos_in_genome = rang_local_ref[dif]
		ntd_in_MCS = MCS[dif]
		try:
			snp = chromosome[pos_in_genome]
			name_snp = snp[2]
			snps_set.add(name_snp)

			snp_ntds = snp[1].split(',')
			snp_ntds_aux = []

			if strand == '-':
				for ntd in snp_ntds:
					if len(ntd) == 1 and ntd != 'N':
						snp_ntds_aux.append(reverse_translation(ntd))
			else:
				snp_ntds_aux = [ snp_aux for snp_aux in snp_ntds if len(snp_aux) == 1 and snp_aux != 'N']

			if ntd_in_MCS in snp_ntds_aux:
				from_ = snp[0]
				if strand == '-':
					from_ = reverse_translation(snp[0])
				info_snp_to_add = [from_, ntd_in_MCS, name_snp, pos_in_genome, dif]
				list_seq_reference_local[dif] = ntd_in_MCS
				info_snps.append(info_snp_to_add)
			elif not var:
				return MCS_perfect_alignments, MCS_variant_alignments
		except KeyError:
			if not var:
				return MCS_perfect_alignments, MCS_variant_alignments
			else:
				pass

	#if var:
	new_sequence = "".join(list_seq_reference_local)
	local_translation_peptide_aux = translation_seq(chr, new_sequence)
	differences_pep = [peptide[i]+':'+str(i) for i in range(len(peptide)) if peptide[i]!= local_translation_peptide[i]]
	#else:
	#	differences_pep = []
	
	if len(info_snps) == len(differences_ntds):
		MCS_perfect_alignments = [MCS, [local_translation_peptide, differences_pep, info_snps, differences_ntds]]
	elif var and local_translation_peptide_aux == peptide:
		MCS_perfect_alignments = [MCS, [local_translation_peptide, differences_pep, info_snps, differences_ntds]]
	else:
		#if len(differences_ntds) - len(info_snps) <= 2:
		#MCS_variant_alignments = [new_sequence, [local_translation_peptide, differences_pep, info_snps, differences_ntds]]
		MCS_variant_alignments = []

	return MCS_perfect_alignments, MCS_variant_alignments


def get_sequences_at_position_local(peptide, seq_reference_local, MCS, rang_local_ref, strand, chromosome, chr):

	MCS_perfect_alignments = []
	MCS_variant_alignments = []
	
	local_translation_peptide = translation_seq(chr, seq_reference_local)
	differences_pep = [peptide[i]+':'+str(i) for i in range(len(peptide)) if peptide[i]!= local_translation_peptide[i]]
	list_seq_reference_local = list(seq_reference_local)

	info_snps = []
	snps_set = set()
	if len(differences_pep) > 4:
		return MCS_perfect_alignments, MCS_variant_alignments

	for ind, dif in enumerate(differences_pep):
		aa = dif.split(':')[0]
		dif = int(dif.split(':')[1])
		ntds_to_search = rang_local_ref[dif*3:dif*3+3]
		positions_reference = range(dif*3, dif*3+3)
		
		for index, pos_in_genome in enumerate(ntds_to_search):
			ntd_pos = positions_reference[index]
			codon_list = list(seq_reference_local[dif*3:dif*3+3])
			try:
				snp = chromosome[pos_in_genome]
				name_snp = snp[2]
				snps_set.add(name_snp)
				snp_ntds = snp[1].split(',')
				snp_ntds_aux = []

				if strand == '-':
					for ntd in snp_ntds:
						if len(ntd) == 1 and ntd != 'N':
							snp_ntds_aux.append(reverse_translation(ntd))
				else:
					snp_ntds_aux = [ snp_aux for snp_aux in snp_ntds if len(snp_aux) == 1 and snp_aux != 'N']

				into = False
				
				for ntd in snp_ntds_aux:
					codon_list[index] = ntd
					translation_codon = translation_seq(chr, "".join(codon_list))
					if translation_codon == aa:
						from_ = snp[0]
						if strand == '-':
							from_ = reverse_translation(snp[0])
						info_snp_to_add = [from_, ntd, name_snp, pos_in_genome, ntd_pos]
						list_seq_reference_local[ntd_pos] = ntd
						info_snps.append(info_snp_to_add)
						into = True
						break
				if into:
					break
			except KeyError:
				pass

		if not var and (len(info_snps) != ind + 1):
			return MCS_perfect_alignments, MCS_variant_alignments

	new_sequence = "".join(list_seq_reference_local)
	local_translation_peptide_aux = translation_seq(chr, new_sequence)

	try:
		differences_ntds = [new_sequence[i]+':'+str(i) for i in range(len(new_sequence)) if new_sequence[i]!= seq_reference_local[i]]
	except IndexError:
		super_logger.info( 'Index Error %s %s %s %s %s %s %s %s %s ', peptide, MCS, chr, strand, new_sequence, str(len(new_sequence)), seq_reference_local, str(len(seq_reference_local)), str(rang_local_ref))
		differences_ntds = 100

	if peptide == local_translation_peptide_aux and len(info_snps) == len(differences_ntds):
		MCS_perfect_alignments = [new_sequence, [local_translation_peptide, differences_pep, info_snps, differences_ntds]]
	elif var and local_translation_peptide_aux == peptide :
		MCS_perfect_alignments = [new_sequence, [local_translation_peptide, differences_pep, info_snps, differences_ntds]]
	else:
		#if len(differences_ntds) - len(info_snps) <= 2:
		#MCS_variant_alignments = [new_sequence, [local_translation_peptide, differences_pep, info_snps, differences_ntds]]
		MCS_variant_alignments = []
	return MCS_perfect_alignments, MCS_variant_alignments


def translation_seq(chr, seq):
	if chr=='chrM':
		translation = uf.translateDNA(seq, frame = 'f1', translTable_id='mt')
	else:
		translation = uf.translateDNA(seq, frame = 'f1', translTable_id='default')

	return translation


def get_alignments(sam_file, dbSNP, common, super_logger_aux, var_aux, genome_version, mode, mouse):

	global path_to_db
	global super_logger
	global var
	global genomePathFai
	global genomePath
	global mode_translation

	super_logger = super_logger_aux
	var = var_aux

	if genome_version == 'v26_88': 
		genomePathFai = path_to_lib + 'genome_versions/genome_v26_88/GRCh38.primary_assembly.genome.fa.fai'
		genomePath = path_to_lib + 'genome_versions/genome_v26_88/GRCh38.primary_assembly.genome.fa'
	elif genome_version == 'v33_99':
		genomePathFai = path_to_lib + 'genome_versions/genome_v33_99/GRCh38.primary_assembly.genome.fa.fai'
		genomePath = path_to_lib + 'genome_versions/genome_v33_99/GRCh38.primary_assembly.genome.fa'
	else:
		genomePathFai = path_to_lib + 'genome_versions/genome_v38_104/GRCh38.primary_assembly.genome.fa.fai'
		genomePath = path_to_lib + 'genome_versions/genome_v38_104/GRCh38.primary_assembly.genome.fa'

	if mouse:
		if genome_version == 'M24':
			genomePathFai = path_to_lib + 'genome_versions/genome_mouse_m24/GRCm38.primary_assembly.genome.fa.fai'
			genomePath = path_to_lib + 'genome_versions/genome_mouse_m24/GRCm38.primary_assembly.genome.fa'
		if genome_version == 'M30':
			genomePathFai = path_to_lib + 'genome_versions/genome_mouse_m30/GRCm39.primary_assembly.genome.fa.fai'
			genomePath = path_to_lib + 'genome_versions/genome_mouse_m30/GRCm39.primary_assembly.genome.fa'

	if mode == 'translation':
		mode_translation = True
	else:
		mode_translation = False

	super_logger.info('Using genome version %s. ', genomePath)
	
	if mouse:
		if dbSNP == 'mouse_GRCm38':
			path_to_db = path_to_lib+'/snps/snps_dics_mouse_GRCm38/'
		if dbSNP == 'mouse_GRCm39':
			path_to_db = path_to_lib+'/snps/snps_dics_mouse_GRCm39/'
	elif dbSNP == 0 :
		path_to_db = ''
	elif dbSNP == 149:
		if common:
			path_to_db = path_to_lib+'/snps/snps_dics_149_common/'
		else:
			path_to_db = path_to_lib+'/snps/snps_dics_149/'
	elif dbSNP == 151:
		if common:
			path_to_db = path_to_lib+'/snps/snps_dics_151_common/'
		else:
			path_to_db = path_to_lib+'/snps/snps_dics_151/'
	else:
		if common:
			path_to_db = path_to_lib+'/snps/snps_dics_155_common/'
		else:
			path_to_db = path_to_lib+'/snps/snps_dics_155/'


	super_logger.info('Using dbSNP database %s with COMMON SNPs = %s. Database Path : %s ', str(dbSNP), str(common), str(path_to_db))
	
	exists = os.path.exists(sam_file+'.dic')
	if not exists:
		alignments_by_chromosome_strand = read_sam_file(sam_file)
	else:
		with open(sam_file+'.dic', 'rb') as fp:
			alignments_by_chromosome_strand = pickle.load(fp)
		super_logger.info('Information SAM file already collected !')

	od = collections.OrderedDict(sorted(alignments_by_chromosome_strand.items(), reverse=True))

	del alignments_by_chromosome_strand

	positions_mcs_peptides_perfect_alignment = {}
	positions_mcs_peptides_variants_alignment = {}
	total_peptides_in = set()
	
	keys = list(od.keys())
	values = list(od.values())
	nodes = multiprocessing.cpu_count()
	
	if common or dbSNP == 149 or dbSNP == 0:
		pool = ProcessPool(nodes=nodes)
		results = pool.map(get_alignments_chromosome, keys, values)
	
		for res in results:
			positions_mcs_peptides_perfect_alignment.update(res[0])
			positions_mcs_peptides_variants_alignment.update(res[1])
			total_peptides_in = total_peptides_in.union(res[2])

		keys.clear()
		values.clear()

	else:
		for chr in ['chr1', 'chr2', 'chr3', 'chr4', 'chr5']:
			try:
				index = keys.index(chr)
				res = get_alignments_chromosome(chr, values[index])
				positions_mcs_peptides_perfect_alignment.update(res[0])
				positions_mcs_peptides_variants_alignment.update(res[1])
				total_peptides_in = total_peptides_in.union(res[2])
				del keys[index]
				del values[index]
			except ValueError:
				pass
			
		nodes = 5
		cont = 0
		for i in range(0,len(keys),nodes):
			pool = ProcessPool(nodes=nodes)
			cont += nodes
			results = pool.map(get_alignments_chromosome, keys[i:cont], values[i:cont])

			for res in results:
				positions_mcs_peptides_perfect_alignment.update(res[0])
				positions_mcs_peptides_variants_alignment.update(res[1])
				total_peptides_in = total_peptides_in.union(res[2])

			pool.close()
			pool.join()
			pool.clear()

	return positions_mcs_peptides_perfect_alignment, positions_mcs_peptides_variants_alignment, total_peptides_in

