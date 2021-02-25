import warnings, time, os, itertools
warnings.filterwarnings("ignore")
import pysam, re
from itertools import groupby
from operator import itemgetter
import utils.useful_functions as uf
import pickle, logging
import numpy as np
from pathos.multiprocessing import ProcessPool
import multiprocessing
import collections

__author__ = "Maria Virginia Ruiz Cuevas"

path_to_lib = '/'.join(os.path.abspath(__file__).split('/')[:-3])+'/lib/'

genomePathFai = path_to_lib + 'GRCh38.primary_assembly.genome.fa.fai'
genomePath = path_to_lib + 'GRCh38.primary_assembly.genome.fa'


def set_strand_read(strand):
	number = "{0:b}".format(int(strand))
	if len(number)>= 5:
		if '1' == number[-5]:
			return '-'
		else:
			return '+'
	else:
		return '+'

def get_ranges(cigar, start, lenSeq, strand, chr):

	faFile = pysam.FastaFile(genomePath, genomePathFai)
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
	
	if chr=='chrM':
		translation_peptide = uf.translateDNA(seqReference, frame = 'f1', translTable_id='mt')
	else:
		translation_peptide = uf.translateDNA(seqReference, frame = 'f1', translTable_id='default')

	return realEnd-1, rang, operators, seqReference, translation_peptide


def get_local_reference(start, lenSeq, chr, strand):

	faFile = pysam.FastaFile(genomePath, genomePathFai)
	seqReference = ''
	end = start+lenSeq-1
	rang = range(start, end+1, 1)
	seqReference = faFile.fetch(chr, start-1,end)
	if strand == '-':
		seqReference = uf.reverseComplement(seqReference)
		rang = rang[::-1]

	if chr=='chrM':
		local_translation_peptide = uf.translateDNA(seqReference, frame = 'f1', translTable_id='mt')
	else:
		local_translation_peptide = uf.translateDNA(seqReference, frame = 'f1', translTable_id='default')
	return end, rang, seqReference, local_translation_peptide


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
	aligments_by_chromosome_strand = {}
	MCS_aligments_by_chromosome_strand = {}
	
	
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
				number_hit = int(line.strip().split('NH')[1].split(':')[2].split('\t')[0])
				MCS = splitLine[9]
				lenSeq = len(MCS)
				
				if chr.isdigit() or chr == 'X' or chr == 'Y' or chr == 'M' or chr == 'MT':
					if chr == 'MT':
						chr = 'chrM'
					else:
						chr = 'chr'+chr
				if strand == '-':
					MCS = uf.reverseComplement(MCS)

				#Chr7:113332624

				#if chr == 'chr7' and peptide == 'TTRPALQEL' and readStart == 113332624:#== 'chrX' or chr == 'chrY':
				key = str(readStart)+'|'+cigar+'|'+str(lenSeq)+'|'+peptide+'|'+strand
				key_2 = chr+'|'+str(readStart)+'|'+peptide+'|'+strand
				#print (key, key_2)
				try:
					chromosomes_alignments = aligments_by_chromosome_strand[chr]
					chromosomes_alignments.add(key)
				except KeyError:
					set_alignments = set()
					set_alignments.add(key)
					aligments_by_chromosome_strand[chr] = set_alignments
				try:
					set_mcs = MCS_aligments_by_chromosome_strand[key_2]
					set_mcs.add(MCS)
				except KeyError:
					set_mcs = set()
					set_mcs.add(MCS)
					MCS_aligments_by_chromosome_strand[key_2] = set_mcs

	timeFinal = time.time()
	total = (timeFinal-time0) / 60.0
	total_chromosomes = len(aligments_by_chromosome_strand)
	logging.info('Time reading SAM file : %f min Chromosomes %d ', total, total_chromosomes)
	return aligments_by_chromosome_strand, MCS_aligments_by_chromosome_strand


def get_alignments_chromosome(chr, chromosomes_alignments):
	positions_mcs_peptides_perfect_alignment = {}
	positions_mcs_peptides_variants_alignment = {}
	positions_mcs_peptides_out_alignment = {}
	info_mcs_perfect_alignments = ''
	info_mcs_variant_alignments = ''
	info_mcs_out_alignments = ''
	chromosome = {}
	path_to_lib = '/'.join(os.path.abspath(__file__).split('/')[:-3])+'/lib/'

	try:
		with open(path_to_lib+'/snps_dics/dbsnp149_all_'+chr+'.dic', 'rb') as fp:
			chromosome = pickle.load(fp)
	except IOError:
		chromosome = {}

	for position in chromosomes_alignments: 
		snps_information = []
		
		split_position = position.split('|')
		readStart = int(split_position[0])
		cigar = split_position[1]
		lenSeq = int(split_position[2])
		peptide = split_position[3]
		
		if peptide != 'lero':

			strand  = split_position[4]

			end_local_ref, rang_local_ref, seq_reference_local, local_translation_peptide = get_local_reference(readStart, lenSeq, chr, strand)
			positionsTrans = find_ranges(rang_local_ref)
			range_trans_local = chr+':'+positionsTrans

			realEnd, rang, operators, seqReference, pep_translation_seq_reference = get_ranges(cigar, readStart, lenSeq, strand, chr)
			positionsTrans = find_ranges(rang)
			range_trans = chr+':'+positionsTrans
			
			key = peptide+'_'+range_trans
			key_local = peptide+'_'+range_trans_local

			parfait_alignment_local = key_local in positions_mcs_peptides_perfect_alignment.keys()
			variant_alignment_local = key_local in positions_mcs_peptides_variants_alignment.keys()
			out_alignment_local = key_local in positions_mcs_peptides_out_alignment.keys()

			parfait_alignment = key in positions_mcs_peptides_perfect_alignment.keys()
			variant_alignment = key in positions_mcs_peptides_variants_alignment.keys()
			out_alignment = key in positions_mcs_peptides_out_alignment.keys()

			if parfait_alignment or variant_alignment or out_alignment :
				continue
			if key == key_local:
				if local_translation_peptide == peptide:
					positions_mcs_peptides_perfect_alignment[key_local] = [strand, seq_reference_local, local_translation_peptide, [],0]
				else:
					if 'S' not in operators : #'D' not in operators and 'I' not in operators and 
						perfect_sequences_to_return, variant_sequences_to_return, peptide_with_snps_local_reference, seq_alignment_2, out_sequences_to_return = get_snps_information(chr, chromosome, peptide, pep_translation_seq_reference, seqReference, strand, rang)
						if peptide_with_snps_local_reference == peptide:
							positions_mcs_peptides_perfect_alignment[key] = [strand, seq_alignment_2, peptide_with_snps_local_reference, perfect_sequences_to_return,0]

						if len(variant_sequences_to_return) > 0 and '*' not in peptide_with_snps_local_reference:
							positions_mcs_peptides_variants_alignment[key] = [strand, seq_alignment_2, peptide_with_snps_local_reference, variant_sequences_to_return,0]
								
						if len(out_sequences_to_return) > 0:
							positions_mcs_peptides_out_alignment[key] = [strand, seq_alignment_2, peptide_with_snps_local_reference, out_sequences_to_return,0]
			else:
				if local_translation_peptide == peptide:
					positions_mcs_peptides_perfect_alignment[key_local] = [strand, seq_reference_local, local_translation_peptide, [],0]
				else:
					if not (parfait_alignment_local or variant_alignment_local or out_alignment_local):
						perfect_sequences_to_return, variant_sequences_to_return, peptide_with_snps_local_reference, seq_alignment_2, out_sequences_to_return = get_snps_information(chr, chromosome, peptide, local_translation_peptide, seq_reference_local, strand, rang_local_ref)
						if peptide_with_snps_local_reference == peptide:
							positions_mcs_peptides_perfect_alignment[key_local] = [strand, seq_alignment_2, peptide_with_snps_local_reference, perfect_sequences_to_return,0]
					
				if 'S' not in operators : # 'D' not in operators and 'I' not in operators and 
					perfect_sequences_to_return, variant_sequences_to_return, peptide_with_snps_local_reference, seq_alignment_2, out_sequences_to_return = get_snps_information(chr, chromosome, peptide, pep_translation_seq_reference, seqReference, strand, rang)
					if peptide_with_snps_local_reference == peptide:
						positions_mcs_peptides_perfect_alignment[key] = [strand, seq_alignment_2, peptide_with_snps_local_reference, perfect_sequences_to_return,0]

					if len(variant_sequences_to_return) > 0 and '*' not in peptide_with_snps_local_reference:
						positions_mcs_peptides_variants_alignment[key] = [strand, seq_alignment_2, peptide_with_snps_local_reference, variant_sequences_to_return,0]
							
					if len(out_sequences_to_return) > 0:
						positions_mcs_peptides_out_alignment[key] = [strand, seq_alignment_2, peptide_with_snps_local_reference, out_sequences_to_return,0]
		
	return positions_mcs_peptides_perfect_alignment, positions_mcs_peptides_variants_alignment, positions_mcs_peptides_out_alignment


def get_snps_information(chr, chromosome, peptide_origin, peptide_reference, seqReference, strand, rang):

	seq_alignment_2 = seqReference
	consenssus_snp_involved = []
	
	pack_info_snps = []
	to_combine = []
	differences = [ i for i in range(len(peptide_origin)) if peptide_origin[i]!= peptide_reference[i]]

	perfect_sequences_to_return = []
	variant_sequences_to_return = []
	out_sequences_to_return = []

	dic_differences = {}
	differences_info = []
	list_differences = set()
	differences_visited = set()
	
	for dif in differences:
		list_differences.add(peptide_reference[dif]+':'+str(dif))
		codon_dif = dif*3
		ntd_pos_in_read = [codon_dif, codon_dif+1, codon_dif+2]
		pos_ntd_codons = [rang[ntd_pos_in_read[0]], rang[ntd_pos_in_read[1]], rang[ntd_pos_in_read[2]]]
		ntd_codons = [seqReference[ntd_pos_in_read[0]], seqReference[ntd_pos_in_read[1]], seqReference[ntd_pos_in_read[2]]]
		aux_seq = seqReference[ntd_pos_in_read[0]:ntd_pos_in_read[2]+1]
		to_combine = []
		
		for index, ntd in enumerate(ntd_codons):
			pos_in_genome = pos_ntd_codons[index]
			pos_in_read = ntd_pos_in_read[index]
			to_combine_string = aux_seq[index]
			try:
				snp = chromosome[pos_in_genome]
				name_snp = snp[2]
				snp_ntds = snp[1].split(',')

				for snp_ntd in snp_ntds:
					if strand == '-':
						snp_ntd = uf.reverseComplement(snp_ntd)
					to_combine_string+=snp_ntd
					try:
						dic_dif_aux = dic_differences[dif]
						try:
							dic_aux = dic_dif_aux[index]
							dic_aux[snp_ntd] = [ntd, snp_ntd, name_snp, pos_in_genome, pos_in_read]
						except KeyError:
							dic_dif_aux[index] = {snp_ntd:[ntd, snp_ntd, name_snp, pos_in_genome, pos_in_read]}
					except KeyError:
						dic_differences[dif] = {index: {snp_ntd:[ntd, snp_ntd, name_snp, pos_in_genome, pos_in_read]}}
				
				if len(to_combine_string) > 0:
					to_combine.append(to_combine_string)
			except:
				to_combine.append(aux_seq[index])

		c = list(itertools.product(*to_combine))

		for i, pro in enumerate(c):
			s = list(aux_seq)
			dif_dic = []
			if i > 0:
				first_a = 0
				for a, ntd in enumerate(pro):
					s[a] = ntd
					try:
						d = dic_differences[dif][a][ntd]
						dif_dic.append(d)
					except KeyError:
						pass
				aux_seq_tt = "".join(s)
				trans = uf.translateDNA(aux_seq_tt, frame = 'f1')
				if trans == peptide_origin[dif]:
					seq_alignment_2 = seq_alignment_2[:codon_dif] + aux_seq_tt + seq_alignment_2[codon_dif+3:]
					differences_info.append(dif_dic)
					differences_visited.add(peptide_reference[dif]+':'+str(dif))
					break

	peptide_with_snps_local_reference = uf.translateDNA(seq_alignment_2, frame = 'f1')
	dif_differences_visited = list_differences - differences_visited
	final_differences = [i for dif in differences_info for i in dif]
	if len(differences_info) == len(differences):
		perfect_sequences_to_return = [final_differences] + [list_differences] + [dif_differences_visited]
	elif len(differences) - len(differences_info) <=2 :
		variant_sequences_to_return = [final_differences] + [list_differences] + [dif_differences_visited]
	else:
		out_sequences_to_return = [final_differences] + [list_differences] + [dif_differences_visited]

	return perfect_sequences_to_return, variant_sequences_to_return, peptide_with_snps_local_reference, seq_alignment_2, out_sequences_to_return


def save_output_info(self, alignment_types, mcs_info):
	to_write='Peptide'+'\t'+'Strand'+'\t'+'Alignment'+'\t'+'MCS'+'\t'+'Diff_AA'+'\t'+'SNVs'+'\t'+'Mutations_Non_Annotated'+'\n'

	key = str(readStart)+'|'+cigar+'|'+str(lenSeq)+'|'+peptide+'|'+strand
	key_2 = chr+'|'+str(readStart)+'|'+peptide+'|'+strand
	key_local = peptide+'_'+range_trans_local

	for alignment_type in alignment_types:
		for peptide_info, info_alignment in alignment_type.items():
			alignment = peptide_info.split('_')[1]
			peptide = peptide_info.split('_')[0]
			strand = info_alignment[0]
			MCS = info_alignment[1]
			snvs_write = ''
			dif_aa_write = ''
			mutations_write = ''
			try:
				snvs = info_alignment[3][0]
				dif_aas = info_alignment[3][1]
				mutations_non_annotated = info_alignment[3][2]
				for snv in snvs:
					#T->G|snp:rs12036323|GenPos:239044335|MCSPos:10
					snv = '['+snv[0]+'->'+snv[1]+'|snp:'+snv[2]+'|GenPos:'+str(snv[3])+'|MCSPos:'+str(snv[4])+']'
					snvs_write += snv
				for dif_aa in dif_aas:
					dif_aa_write += '['+dif_aa+']'
				for mutation in mutations_non_annotated:
					mutation_write = '['+mutation+']'
					mutations_write += mutation_write
			except IndexError:
				pass
			to_write+=peptide+'\t'+strand+'\t'+alignment+'\t'+MCS+'\t'+dif_aa_write+'\t'+snvs_write+'\t'+mutations_write+'\n'

	file_to_open = open(self.path_to_save+self.name_exp+name_file, 'w')
	file_to_open.write(to_write)
	file_to_open.close()
	logging.info('Info perfect alignments saved to : %s ', self.path_to_save+self.name_exp+name_file)


def get_alignments(sam_file):

	aligments_by_chromosome_strand, MCS_aligments_by_chromosome_strand = read_sam_file(sam_file)
	
	od = collections.OrderedDict(sorted(aligments_by_chromosome_strand.items(), reverse=True))

	positions_mcs_peptides_perfect_alignment = {}
	positions_mcs_peptides_variants_alignment = {}
	positions_mcs_peptides_out_alignment = {}
	
	nodes = multiprocessing.cpu_count()
	
	keys = od.keys()
	values = od.values()
	pool = ProcessPool(nodes=nodes)
	results = pool.map(get_alignments_chromosome, keys, values)
	for res in results:
		positions_mcs_peptides_perfect_alignment.update(res[0])
		positions_mcs_peptides_variants_alignment.update(res[1])
		positions_mcs_peptides_out_alignment.update(res[2])

	return positions_mcs_peptides_perfect_alignment, positions_mcs_peptides_variants_alignment, positions_mcs_peptides_out_alignment

