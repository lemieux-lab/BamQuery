import warnings
warnings.filterwarnings("ignore")
import time, os
import pysam, re, sys
from itertools import groupby
from operator import itemgetter
import utils.useful_functions as uf
import pickle
from pathos.multiprocessing import ProcessPool
import collections
import pandas as pd

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


def get_ranges(cigar, start, strand, seqReference):
    
    rang = []
    seqReference = seqReference.upper()
	# https://pysam.readthedocs.io/en/latest/api.html?highlight=read.cigartuples#pysam.AlignedSegment.cigartuples
    for tuple_ in cigar:
        operation = tuple_[0]
        length = tuple_[1]

        if (4 == operation): 
            end = length
            rang.extend([0]*length)

        elif (0 == operation) or (7 == operation) or (8 == operation) :
            end = start + length
            rang.extend(range(start,end))
            start += length

        elif (3 == operation) :
            end = start + length
            start = end

    if strand == '-':
        seqReference = uf.reverseComplement(seqReference)
        rang = rang[::-1]

    if 'N' in seqReference:
        seqReference = seqReference.replace("N", "C")
    if 'n' in seqReference:
        seqReference = seqReference.replace("n", "C")
    
    return rang, seqReference

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


def find_ranges(chr, iterable, strand):
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
	return chr+':'+toReturn


def read_sam_file(sam_file):
	time0 = time.time()
	alignments_by_chromosome_strand = {}
	
	samfile = pysam.AlignmentFile(sam_file, "r")
	for read in samfile:
		queryname = read.query_name
		cigar = read.cigarstring

		if 'I' in cigar or 'D' in cigar:
			continue
		cigar = read.cigartuples

		chr = read.reference_name
		readStart = int(read.reference_start) + 1
		seqReference = read.get_reference_sequence()

		if read.is_reverse:
			strand = "-"
		else:
			strand = "+"
			
		peptide = queryname
		MCS = read.query_sequence
		
		if chr.isdigit() or chr == 'X' or chr == 'Y' or chr == 'M' or chr == 'MT':
			if chr == 'MT':
				chr = 'chrM'
			else:
				chr = 'chr'+chr

		if strand == '-':
			MCS = uf.reverseComplement(MCS)

		key = [readStart, cigar, MCS, peptide, strand, seqReference]
		
		try:
			chromosomes_alignments = alignments_by_chromosome_strand[chr]
			chromosomes_alignments.append(key)
		except KeyError:
			set_alignments = []
			set_alignments.append(key)
			alignments_by_chromosome_strand[chr] = set_alignments

	samfile.close()		
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
		readStart = int(position[0])
		cigar = position[1]
		MCS = position[2]
		peptide = position[3]
		strand  = position[4]
		seqReference = position[5]
		lenSeq = len(MCS)
		
		end = readStart+lenSeq-1
		range_trans_local = find_ranges(chr, range(readStart, end+1), strand)
		key_local = peptide+'_'+range_trans_local
		
		rang, seq_reference_align = get_ranges(cigar, readStart, strand, seqReference)
		range_trans = find_ranges(chr, rang, strand)
		if range_trans == chr:
			continue
		key_pos = peptide+'_'+range_trans
		
		MCS_perfect_alignments_align, MCS_variant_alignments_align = get_sequences_at_position(peptide, seq_reference_align, MCS, rang, strand, chromosome, chr)
		
		if len(MCS_perfect_alignments_align) > 0:
			peptides_in.add(peptide)

			MCS_in =  MCS_perfect_alignments_align[0]
			key_p = key_pos+'_'+MCS_in

			local_translation_peptide = MCS_perfect_alignments_align[1][0]
			differences_pep = MCS_perfect_alignments_align[1][1]
			info_snps = MCS_perfect_alignments_align[1][2]
			differences_ntds = MCS_perfect_alignments_align[1][3]

			positions_mcs_peptides_perfect_alignment[key_p] = [strand, local_translation_peptide, differences_pep, info_snps, differences_ntds, []]
			

		if len(MCS_variant_alignments_align) > 0:
			peptides_in.add(peptide)
			MCS_in_var =  MCS_variant_alignments_align[0]
			key_v = key_pos+'_'+MCS_in_var

			local_translation_peptide = MCS_variant_alignments_align[1][0]
			differences_pep = MCS_variant_alignments_align[1][1]
			info_snps = MCS_variant_alignments_align[1][2]
			differences_ntds = MCS_variant_alignments_align[1][3]

			positions_mcs_peptides_perfect_alignment[key_v] = [strand, local_translation_peptide, differences_pep, info_snps, differences_ntds, []]
			

		if key_local not in local_visited and key_pos != key_local:
			local_visited.add(key_local)
			rang_local_ref, seq_reference_local = get_local_reference(readStart, lenSeq, chr, strand, faFile)
			MCS_perfect_alignments_local, MCS_variant_alignments_local = get_sequences_at_position_local(peptide, seq_reference_local, MCS, list(rang_local_ref), strand, chromosome, chr)
			
			if len(MCS_perfect_alignments_local) > 0:
				peptides_in.add(peptide)
			
				MCS_in =  MCS_perfect_alignments_local[0]
				key_local_p = key_local+'_'+MCS_in
				
				local_translation_peptide = MCS_perfect_alignments_local[1][0]
				differences_pep = MCS_perfect_alignments_local[1][1]
				info_snps = MCS_perfect_alignments_local[1][2]
				differences_ntds = MCS_perfect_alignments_local[1][3]

				if key_local_p not in positions_mcs_peptides_perfect_alignment.keys():
					positions_mcs_peptides_perfect_alignment[key_local_p] = [strand, local_translation_peptide, differences_pep, info_snps, differences_ntds, []]
				
			if len(MCS_variant_alignments_local) > 0:
				peptides_in.add(peptide)
				
				MCS_in_var =  MCS_variant_alignments_local[0]
				key_local_v = key_local+'_'+MCS_in_var

				local_translation_peptide = MCS_variant_alignments_local[1][0]
				differences_pep = MCS_variant_alignments_local[1][1]
				info_snps = MCS_variant_alignments_local[1][2]
				differences_ntds = MCS_variant_alignments_local[1][3]

				if key_local_v not in positions_mcs_peptides_perfect_alignment.keys():
					positions_mcs_peptides_perfect_alignment[key_local_v] = [strand, local_translation_peptide, differences_pep, info_snps, differences_ntds, []]

	chromosome = {}
	return positions_mcs_peptides_perfect_alignment, peptides_in


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
	
	if MCS == seq_reference_local:
		MCS_perfect_alignments = [MCS, [peptide, [], [], []]]
		return MCS_perfect_alignments, MCS_variant_alignments

	local_translation_peptide = translation_seq(chr, seq_reference_local)
	differences_ntds = []
	for i in range(len(seq_reference_local)) :
		if seq_reference_local[i]!= MCS[i]:
			differences_ntds.append(seq_reference_local[i]+':'+str(i))
			
	list_seq_reference_local = list(seq_reference_local)
	
	info_snps = []
	snps_set = set()
	
	if len(differences_ntds) > 4:
		return MCS_perfect_alignments, MCS_variant_alignments

	difs_ntd = set()
	for dif in differences_ntds:
		dif = int(dif.split(':')[1])
		difs_ntd.add(dif)
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
	
	new_sequence = "".join(list_seq_reference_local)
	local_translation_peptide_aux = translation_seq(chr, new_sequence)
	differences_pep = []
	difs_aa = []
	for i in range(len(peptide)) :
		if peptide[i]!= local_translation_peptide[i]:
			differences_pep.append(peptide[i]+':'+str(i))
			difs_aa.append(i)
	
	if MCS == new_sequence:
		MCS_perfect_alignments = [MCS, [local_translation_peptide, differences_pep, info_snps, differences_ntds]]
	
	if var and local_translation_peptide_aux == peptide:
		if len(info_snps) == 0 :
			MCS_perfect_alignments = [new_sequence, [local_translation_peptide, [], [], []]]
		else:
			MCS_perfect_alignments = [new_sequence, [local_translation_peptide, differences_pep, info_snps, differences_ntds]]
	
		if (len(differences_ntds) + len(info_snps)) <= 3:
			MCS_variant_alignments = [MCS, [local_translation_peptide, differences_pep, info_snps, differences_ntds]]

	if len(MCS_variant_alignments) == 0 and var and len(difs_aa) == 1  :
		rang = set(range(difs_aa[0]*3, difs_aa[0]*3+3))
		if len(rang.intersection(difs_ntd)) == len(difs_ntd):
			MCS_variant_alignments = [MCS, [local_translation_peptide, differences_pep, [], differences_ntds]]
	
	return MCS_perfect_alignments, MCS_variant_alignments


def get_sequences_at_position_local(peptide, seq_reference_local, MCS, rang_local_ref, strand, chromosome, chr):

	MCS_perfect_alignments = []
	MCS_variant_alignments = []
	
	local_translation_peptide = translation_seq(chr, seq_reference_local)
	differences_pep = [peptide[i]+':'+str(i) for i in range(len(peptide)) if peptide[i]!= local_translation_peptide[i]]
	list_seq_reference_local = list(seq_reference_local)

	info_snps = []
	snps_set = set()
	if len(differences_pep) > 2:
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
	differences_ntds = [new_sequence[i]+':'+str(i) for i in range(len(new_sequence)) if new_sequence[i]!= seq_reference_local[i]]
	
	if peptide == local_translation_peptide_aux and len(info_snps) == len(differences_ntds):
		MCS_perfect_alignments = [new_sequence, [local_translation_peptide, differences_pep, info_snps, differences_ntds]]
		
	return MCS_perfect_alignments, MCS_variant_alignments

def translation_seq(chr, seq):
	if chr=='chrM':
		translation = uf.translateDNA(seq, frame = 'f1', translTable_id='mt')
	else:
		translation = uf.translateDNA(seq, frame = 'f1', translTable_id='default')

	return translation


def get_alignments(sam_file, dbSNP, common, super_logger_aux, var_aux, genome_version, mode, mouse, threads):

	global path_to_db
	global super_logger
	global var
	global genomePathFai
	global genomePath
	global mode_translation
	global splice_junctions_annotated

	super_logger = super_logger_aux
	var = var_aux
	
	if genome_version == 'v26_88': 
		genomePathFai = path_to_lib + 'genome_versions/genome_v26_88/GRCh38.primary_assembly.genome.fa.fai'
		genomePath = path_to_lib + 'genome_versions/genome_v26_88/GRCh38.primary_assembly.genome.fa'
		splice_junctions = path_to_lib + 'genome_versions/genome_v26_88/Index_STAR_2.7.9a/sjdbList.fromGTF.out.tab'
	elif genome_version == 'v33_99':
		genomePathFai = path_to_lib + 'genome_versions/genome_v33_99/GRCh38.primary_assembly.genome.fa.fai'
		genomePath = path_to_lib + 'genome_versions/genome_v33_99/GRCh38.primary_assembly.genome.fa'
		splice_junctions = path_to_lib + 'genome_versions/genome_v33_99/Index_STAR_2.7.9a/sjdbList.fromGTF.out.tab'
	else:
		genomePathFai = path_to_lib + 'genome_versions/genome_v38_104/GRCh38.primary_assembly.genome.fa.fai'
		genomePath = path_to_lib + 'genome_versions/genome_v38_104/GRCh38.primary_assembly.genome.fa'
		splice_junctions = path_to_lib + 'genome_versions/genome_v38_104/Index_STAR_2.7.9a/sjdbList.fromGTF.out.tab'

	if mouse:
		if genome_version == 'M24':
			genomePathFai = path_to_lib + 'genome_versions/genome_mouse_m24/GRCm38.primary_assembly.genome.fa.fai'
			genomePath = path_to_lib + 'genome_versions/genome_mouse_m24/GRCm38.primary_assembly.genome.fa'
			splice_junctions = path_to_lib + 'genome_versions/genome_mouse_m24/Index_STAR_2.7.9a/sjdbList.fromGTF.out.tab'
		if genome_version == 'M30':
			genomePathFai = path_to_lib + 'genome_versions/genome_mouse_m30/GRCm39.primary_assembly.genome.fa.fai'
			genomePath = path_to_lib + 'genome_versions/genome_mouse_m30/GRCm39.primary_assembly.genome.fa'
			splice_junctions = path_to_lib + 'genome_versions/genome_mouse_m30/Index_STAR_2.7.9a/sjdbList.fromGTF.out.tab'

	splice_junctions_annotated = pd.read_csv(splice_junctions, header=None, sep='\t')

	if mode == 'translation':
		mode_translation = True
	else:
		mode_translation = False

	super_logger.info('Using genome version %s. ', genomePath)
	
	if mouse:
		if dbSNP == 'mouse_GRCm38':
			path_to_db = path_to_lib+'/snps/snps_dics_mouse_GRCm38/'
			exists_db = os.path.exists(path_to_db) 
			if not exists_db:
				path_to_db = ''
				super_logger.info('dbSNPs for the M24 mouse genome is not found in the lib folder. If you want to use snps, follow the instructions at https://bamquery.iric.ca/documentation/installation.html to download snps_mouse_GRCm38.')
		if dbSNP == 'mouse_GRCm39':
			path_to_db = path_to_lib+'/snps/snps_dics_mouse_GRCm39/'
			exists_db = os.path.exists(path_to_db) 
			if not exists_db:
				path_to_db = ''
				super_logger.info('dbSNPs for the M30 mouse genome is not found in the lib folder. If you want to use snps, follow the instructions at https://bamquery.iric.ca/documentation/installation.html to download snps_mouse_GRCm39.')
		
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
	total_peptides_in = set()
	
	keys = list(od.keys())
	values = list(od.values())
	
	if common or dbSNP == 149 or dbSNP == 0:
		pool = ProcessPool(nodes=threads)
		results = pool.map(get_alignments_chromosome, keys, values)
	
		for res in results:
			positions_mcs_peptides_perfect_alignment.update(res[0])
			total_peptides_in = total_peptides_in.union(res[1])

		keys.clear()
		values.clear()

	else:
		for chr in ['chr1', 'chr2', 'chr3', 'chr4', 'chr5']:
			try:
				index = keys.index(chr)
				res = get_alignments_chromosome(chr, values[index])
				positions_mcs_peptides_perfect_alignment.update(res[0])
				total_peptides_in = total_peptides_in.union(res[1])
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
				total_peptides_in = total_peptides_in.union(res[1])

			pool.close()
			pool.join()
			pool.clear()
			
	return positions_mcs_peptides_perfect_alignment, total_peptides_in

