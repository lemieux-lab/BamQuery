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

	bam_handle = pysam.AlignmentFile(bam_file, mode = 'rb')

	ires = map(
		_get_counts_sample,
		repeat(name_sample),
		repeat(library),
		repeat(sens),
		repeat(bam_handle),
		peptide_alignment,
		sequences,
		repeat(references),
		repeat(overlap)
	)

	res = [list(x) for x in ires]

	bam_handle.close()

	return res

def _get_counts_sample(
		name_sample,
		library,
		sens,
		bam_handle,
		peptide_alignment,
		sequences,
		references,
		overlap
	):

	count = 0
	peptide = peptide_alignment.split('_')[0]
	alignment = peptide_alignment.split('_')[1]
	strand = peptide_alignment.split('_')[2]

	yield [peptide, alignment, name_sample, strand]

	chr = alignment.split(':')[0]
	try:
		chr = references[chr]
	except:
		for sequence in sequences:
			yield [0, sequence]
		return

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

	counts = get_depth_with_view(
		region_to_query,
		bam_handle,
		library,
		sens,
		strand,
		sequences,
		pos_set,
		splice_pos,
		overlap
	)

	for index, sequence in enumerate(sequences):
		try:
			count = counts[index]
		except TypeError:
			count = -1
		yield [count, sequence]


def set_strand_read(strand):
	number = "{0:b}".format(int(strand))
	if len(number)>= 5:
		if '1' == number[-5]:
			return '-'
		else:
			return '+'
	else:
		return '+'

def get_ranges_aux(cigar, start, len_seq):

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

def get_ranges(cigar, start, len_seq):

	rang = []
	splices_sites = set()
	#rx = re.findall('(\d+)([MISDNX=])?', cigar)
	# M == 0
	# I == 1
	# D == 2
	# N == 3
	# S == 4
	# = == 7
	# X == 8
	indels = []

	for operation, length in cigar:
		if operation == 4:
			end = length
			rang.extend([0]*length)

		elif operation == 0 or operation == 7 or operation == 8:
			end = start + length
			rang.extend(range(start,end))
			start += length

		elif operation == 3 or operation == 2:
			splices_sites.add(start)
			start = start + length
			splices_sites.add(start)

		elif operation == 1:
			indels.extend([len(rang)+i for i in range(0, length)])
			rang.extend([0]*length)

	return rang, splices_sites, indels

def remove_at(i, s):
	return s[:i] + s[i+1:]

def define_read_in(reads, pos_set, splice_pos, overlap_in):
	for read in reads:
		name = read.query_name
		strand = '+' if read.is_forward else '-'
		start = read.reference_start + 1
		seq = read.query_sequence
		cigar = tuple(read.cigartuples)

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
			seq_overlap = seq[index_ini:index_fin]
			if overlap_in and percentage_overlap >= 0.6 :
				percentage_overlap = 1.0

			if percentage_overlap == 1.0:
				yield name, strand, seq_overlap, percentage_overlap

def get_depth_with_view(
		region_to_query,
		bam,
		library,
		sens,
		strand,
		sequences,
		pos_set,
		splice_pos,
		overlap
	):

	contReads_to_return = []

	library = library.lower()
	sens = sens.lower()

	def set_count(info_read, cs):

		name, strand, seq_overlap, percentage_overlap = info_read

		if seq_overlap in cs and name not in set_names_reads:
			set_names_reads.add(name)
			return percentage_overlap

	try:
		reads = bam.fetch(region=region_to_query)
	except pysam.utils.SamtoolsError:
		return -1

	reads = (read for read in reads if not read.is_unmapped)
	# equivalent to samtools view flag '-F0X4'

	if library == 'unstranded' or sens == 'unstranded':
		reads = (read for read in reads if not read.is_secondary)
		# equivalent to samtools view flag '-F0X100'

		reads_overlaping_area = list(define_read_in(reads, pos_set, splice_pos, overlap))
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
			reads = (read for read in reads if not read.is_secondary and read.is_forward)
			# equivalent to samtools view flag '-F0X110'

		elif ((strand == '-' and sens == 'forward') or (strand == '+' and sens == 'reverse')):
			reads = (read for read in reads if not read.is_secondary and read.is_reverse)
			# equivalent to samtools view flags '-F0X100' and '-f0X10'

		reads_overlaping_area = list(define_read_in(reads, pos_set, splice_pos, overlap))
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
		reads = [read for read in reads if not read.is_secondary]
		# equivalent to samtools view flag '-F0X100'

		if ((strand == '+' and sens == 'forward') or (strand == '-' and sens == 'reverse')):
			reads_1 = (read for read in reads if read.is_read1 and read.mate_is_reverse)
			# equivalent to samtools view flag '-f0X60'
			reads_2 = (read for read in reads if read.is_read2 and read.is_reverse)
			# equivalent to samtools view flag '-f0X90'

		elif ((strand == '-' and sens == 'forward') or (strand == '+' and sens == 'reverse')):
			reads_1 = (read for read in reads if read.is_read1 and read.is_reverse)
			# equivalent to samtools view flag '-f0X50'
			reads_2 = (read for read in reads if read.is_read2 and read.mate_is_reverse)
			# equivalent to samtools view flag '-f0XA0'

		reads_overlaping_area_count_1 = list(define_read_in(reads_1, pos_set, splice_pos, overlap))
		reads_overlaping_area_count_1.sort(key=itemgetter(-1), reverse=True)

		reads_overlaping_area_count_2 = list(define_read_in(reads_2, pos_set, splice_pos, overlap))
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
