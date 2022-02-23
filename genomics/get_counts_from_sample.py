import utils.useful_functions as uf
import pysam


__author__ = "Maria Virginia Ruiz Cuevas"
__email__ = "maria.virginia.ruiz.cuevas@umontreal.ca"

def get_counts_sample(bam, peptide_alignment, sequences):
		
	name_sample = bam[0]
	bam_file = bam[1][0]
	library = bam[1][1]
	sens = bam[1][2]

	count = 0
	peptide = peptide_alignment.split('_')[0]
	alignment = peptide_alignment.split('_')[1]
	strand = peptide_alignment.split('_')[2]

	chr = alignment.split(':')[0]
	
	region_to_query = chr+':'+alignment.split(':')[1].split('-')[0]+'-'+alignment.split(':')[1].split('-')[-1]

	counts = get_depth_with_view(region_to_query, bam_file, library, sens, strand, sequences)

	to_return = [[peptide, alignment, name_sample, strand]]
	
	for index, sequence in enumerate(sequences):
		try:
			count = counts[index]
		except TypeError:
			count = -1
		to_return.append([count, sequence])
	 
	return to_return


def get_depth_with_view(region_to_query, bam_file, library, sens, strand, sequences):

	contReads_to_return = []

	library = library.lower()
	sens = sens.lower()
	
	if library == 'unstranded' or sens == 'unstranded':
		try:
			count_1 = pysam.view("-F0X100", bam_file, region_to_query)
		except pysam.utils.SamtoolsError: 
			return -1

		for seq in sequences:
			contReads = 0
			rcmcs = uf.reverseComplement(seq)
			contReads += count_1.count(seq)
			contReads += count_1.count(rcmcs)
			contReads_to_return.append(contReads)

	elif library == 'single-end':

		if ((strand == '+' and sens == 'forward') or (strand == '-' and sens == 'reverse')):
			try:
				count_1 = pysam.view("-F0X110", bam_file, region_to_query)
			except pysam.utils.SamtoolsError: 
				return -1
		elif ((strand == '-' and sens == 'forward') or (strand == '+' and sens == 'reverse')):
			try:
				count_1 = pysam.view("-F0X100", "-f0X10", bam_file, region_to_query)
			except pysam.utils.SamtoolsError: 
				return -1

		for seq in sequences:
			contReads = 0
			contReads += count_1.count(seq)
			contReads_to_return.append(contReads)
			
	elif library == 'pair-end':
		if ((strand == '+' and sens == 'forward') or (strand == '-' and sens == 'reverse')):
			try:
				count_1 = pysam.view("-F0X100", "-f0X60", bam_file, region_to_query)
				count_2 = pysam.view("-F0X100", "-f0X90", bam_file, region_to_query)
			except pysam.utils.SamtoolsError: 
				return -1

		elif ((strand == '-' and sens == 'forward') or(strand == '+' and sens == 'reverse')):
			try:
				count_1 = pysam.view("-F0X100", "-f0X50", bam_file, region_to_query)
				count_2 = pysam.view("-F0X100", "-f0XA0", bam_file, region_to_query)
			except pysam.utils.SamtoolsError: 
				return -1

		count_1_split = count_1.split('\n')
		count_2_split = count_2.split('\n')
		
		for seq in sequences:
			contReads = 0
			reads_name = set()
			rcmcs = uf.reverseComplement(seq)

			for read in count_1_split:
				if seq in read :
					name = read.split('\t')[0]
					reads_name.add(name)

			for read in count_2_split:
				if rcmcs in read :
					name = read.split('\t')[0]
					reads_name.add(name)

			contReads = len(reads_name)
			contReads_to_return.append(contReads)

		# reading htslib https://github.com/DecodeGenetics/graphtyper/issues/57
	return contReads_to_return
	