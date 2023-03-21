import utils.useful_functions as uf
import pysam, time
import re
from itertools import repeat
from operator import itemgetter

__author__ = "Maria Virginia Ruiz Cuevas"
__email__ = "maria.virginia.ruiz.cuevas@umontreal.ca"


def get_counts_sample(bam, peptide_alignment, sequences, references, overlap):
		
	name_sample = bam[0]
	bam_file = bam[1][0]
	library = bam[1][1]
	sens = bam[1][2]

	count = 0
	peptide = peptide_alignment.split('_')[0]
	alignment = peptide_alignment.split('_')[1]
	strand = peptide_alignment.split('_')[2]
	to_return = [[peptide, alignment, name_sample, strand]]

	chr = alignment.split(':')[0]
	try:
		chr = references[chr]
	except:
		for sequence in sequences:
			to_return.append([0, sequence])
		return to_return

	region_to_query = chr+':'+alignment.split(':')[1].split('-')[0]+'-'+alignment.split(':')[1].split('-')[-1]

	pos = alignment.split(':')[1].split('|')
	pos_set = []
	splice_pos = set()
	for i, chunck in enumerate(pos):
		ini = int(chunck.split('-')[0])
		fini = int(chunck.split('-')[1])
		aux = set(range(ini, fini+1))
		pos_set.extend(aux)
	
		if i != 0:
			splice_pos.add(ini)
		if len(pos) > 1 and i != len(pos) -1:
			splice_pos.add(fini)

	counts = get_depth_with_view(region_to_query, bam_file, library, sens, strand, sequences, pos_set, splice_pos, overlap)
	
	for index, sequence in enumerate(sequences):
		try:
			count = counts[index]
		except TypeError:
			count = -1
		to_return.append([count, sequence])
	return to_return


def set_strand_read(strand):
	number = "{0:b}".format(int(strand))
	if len(number)>= 5:
		if '1' == number[-5]:
			return '-'
		else:
			return '+'
	else:
		return '+'

def get_ranges(cigar, start, len_seq):
	
	rang = [0]*len_seq
	splices_sites = set()
	rx = re.findall('(\d+)([MISDNX=])?', cigar)
	indels = []
	lastIndex = 0

	for index in rx:
		operation = index[1]
		length = int(index[0])
		
		if ('S' in operation):
			end = length
			lastIndex += end
			
		elif ('M' in index) or ('=' in index) or ('X' in index) :
			end = length

			for i in range(0,end):
				rang[lastIndex+i] = start
				start = start + 1

			lastIndex = lastIndex + end
			
		elif ('N' in index) or ('D' in index):
			splices_sites.add(start)
			start = start + length
			splices_sites.add(start)

		elif ('I' in index):
			indels.append(len(rang)+1)

	return rang, splices_sites, indels

def remove_at(i, s):
	return s[:i] + s[i+1:]

def define_read_in(read, pos_set, splice_pos, overlap_in):

	if len(read) > 0:
		split_read = read.split('\t')
		name = split_read[0]
		strand = split_read[1]
		chr = split_read[2]
		start = int(split_read[3])
		seq = split_read[9]
		cigar = split_read[5]
		strand = set_strand_read(strand)
		rang_, splices_sites, indels = get_ranges(cigar, start, len(seq))
		splices_sites = splices_sites - splice_pos
		overlap = set(rang_).intersection(pos_set)
		percentage_overlap = len(overlap)/len(pos_set)
		
		if len(indels) > 0:
			for indel in indels:
				seq = remove_at(indel, seq)
		
		if len(splices_sites.intersection(sorted(pos_set)[1:-1])) == 0 and percentage_overlap >= 0.6:
			index_ini = rang_.index(min(overlap))
			index_fin = rang_.index(max(overlap)) + 1
			seq_overlap = seq[index_ini :index_fin]
			if overlap_in:
				if percentage_overlap >= 0.6 :
					return name, strand, seq_overlap, percentage_overlap
			else:
				if percentage_overlap == 1.0 :
					return name, strand, seq_overlap, percentage_overlap
		

def get_depth_with_view(region_to_query, bam_file, library, sens, strand, sequences, pos_set, splice_pos, overlap):

	contReads_to_return = []

	library = library.lower()
	sens = sens.lower()
	
	def set_count(info_read, cs):

		name, strand, seq_overlap, percentage_overlap = info_read
		
		if seq_overlap in cs and name not in set_names_reads:
			set_names_reads.add(name)
			return percentage_overlap 

	if library == 'unstranded' or sens == 'unstranded':
		try:
			count_1 = pysam.view("-F0X100", bam_file, region_to_query).split('\n')
		except pysam.utils.SamtoolsError: 
			return -1

		reads_overlaping_area = list(map(define_read_in, count_1, repeat(pos_set), repeat(splice_pos), repeat(overlap)))
		reads_overlaping_area = list(filter(None, reads_overlaping_area))
		reads_overlaping_area.sort(key=itemgetter(-1), reverse=True)

		for seq in sequences:
			contReads = 0
			set_names_reads = set()
			rcmcs = uf.reverseComplement(seq)
			
			sum_overlap_seq = list(map(set_count, reads_overlaping_area, repeat(seq)))
			sum_overlap_seq = sum(list(filter(None, sum_overlap_seq)))

			sum_overlap_rcmcs= list(map(set_count, reads_overlaping_area, repeat(rcmcs)))
			sum_overlap_rcmcs = sum(list(filter(None, sum_overlap_rcmcs)))
			
			contReads = sum_overlap_seq + sum_overlap_rcmcs
			
			contReads_to_return.append(contReads)

	elif library == 'single-end':

		if ((strand == '+' and sens == 'forward') or (strand == '-' and sens == 'reverse')):
			try:
				count_1 = pysam.view("-F0X110", bam_file, region_to_query).split('\n')
			except pysam.utils.SamtoolsError: 
				return -1
		elif ((strand == '-' and sens == 'forward') or (strand == '+' and sens == 'reverse')):
			try:
				count_1 = pysam.view("-F0X100", "-f0X10", bam_file, region_to_query).split('\n')
			except pysam.utils.SamtoolsError: 
				return -1

		reads_overlaping_area = list(map(define_read_in, count_1, repeat(pos_set), repeat(splice_pos), repeat(overlap)))
		reads_overlaping_area = list(filter(None, reads_overlaping_area))
		reads_overlaping_area.sort(key=itemgetter(-1), reverse=True)

		for seq in sequences:
			if strand == '-':
				seq = uf.reverseComplement(seq)

			contReads = 0
			set_names_reads = set()
			sum_overlap_seq = list(map(set_count, reads_overlaping_area, repeat(seq)))
			contReads = sum(list(filter(None, sum_overlap_seq)))

			contReads_to_return.append(contReads)
			
	elif library == 'pair-end':
		if ((strand == '+' and sens == 'forward') or (strand == '-' and sens == 'reverse')):
			try:
				count_1 = pysam.view("-F0X100", "-f0X60", bam_file, region_to_query).split('\n')
				count_2 = pysam.view("-F0X100", "-f0X90", bam_file, region_to_query).split('\n')
			except pysam.utils.SamtoolsError: 
				return -1

		elif ((strand == '-' and sens == 'forward') or (strand == '+' and sens == 'reverse')):
			try:
				count_1 = pysam.view("-F0X100", "-f0X50", bam_file, region_to_query).split('\n')
				count_2 = pysam.view("-F0X100", "-f0XA0", bam_file, region_to_query).split('\n')
			except pysam.utils.SamtoolsError: 
				return -1

		reads_overlaping_area_count_1 = list(map(define_read_in, count_1, repeat(pos_set), repeat(splice_pos), repeat(overlap)))
		reads_overlaping_area_count_1 = list(filter(None, reads_overlaping_area_count_1))
		reads_overlaping_area_count_1.sort(key=itemgetter(-1), reverse=True)

		reads_overlaping_area_count_2 = list(map(define_read_in, count_2, repeat(pos_set), repeat(splice_pos), repeat(overlap)))
		reads_overlaping_area_count_2 = list(filter(None, reads_overlaping_area_count_2))
		reads_overlaping_area_count_2.sort(key=itemgetter(-1), reverse=True)

		for seq in sequences:
			contReads = 0
			set_names_reads = set()
			rcmcs_aux = uf.reverseComplement(seq)

			if strand == '-':
				seq = rcmcs_aux
				rcmcs = rcmcs_aux
			else:
				rcmcs = seq

			sum_overlap_seq = list(map(set_count, reads_overlaping_area_count_1, repeat(seq)))
			sum_overlap_seq = sum(list(filter(None, sum_overlap_seq)))

			sum_overlap_rcmcs= list(map(set_count, reads_overlaping_area_count_2, repeat(rcmcs)))
			sum_overlap_rcmcs = sum(list(filter(None, sum_overlap_rcmcs)))

			contReads = sum_overlap_seq + sum_overlap_rcmcs
			contReads_to_return.append(contReads)

	# reading htslib https://github.com/DecodeGenetics/graphtyper/issues/57
	return contReads_to_return

	