import warnings, time, logging, os
warnings.filterwarnings("ignore")
import concurrent.futures
import multiprocessing

NUM_WORKERS =  multiprocessing.cpu_count()

CODON_TABLE = {
	'A': ('GCT', 'GCC', 'GCA', 'GCG'),
	'C': ('TGT', 'TGC'),
	'D': ('GAT', 'GAC'),
	'E': ('GAA', 'GAG'),
	'F': ('TTT', 'TTC'),
	'G': ('GGT', 'GGC', 'GGA', 'GGG'),
	'I': ('ATT', 'ATC', 'ATA'),
	'H': ('CAT', 'CAC'),
	'K': ('AAA', 'AAG'),
	'L': ('TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'),
	'M': ('ATG',),
	'N': ('AAT', 'AAC'),
	'P': ('CCT', 'CCC', 'CCA', 'CCG'),
	'Q': ('CAA', 'CAG'),
	'R': ('CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'),
	'S': ('TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'),
	'T': ('ACT', 'ACC', 'ACA', 'ACG'),
	'V': ('GTT', 'GTC', 'GTA', 'GTG'),
	'W': ('TGG',),
	'Y': ('TAT', 'TAC'),
	'*': ('TAA', 'TAG', 'TGA'),
}

__author__ = "Maria Virginia Ruiz Cuevas"

class ReverseTranslation:

	def __init__(self):
		pass

	def reverse_translation(self, peptide_mode, CS_mode, path_to_output_folder, name_exp):
		
		exist = os.path.exists(path_to_output_folder+'genome_alignments/'+name_exp+'.fastq')
		if not exist:
			t_0 = time.time()
			to_write_reverse_translation = ''
			total_reads_to_align = 0
			
			for peptide in peptide_mode:
				with concurrent.futures.ThreadPoolExecutor(max_workers=NUM_WORKERS) as executor:
					future = executor.submit(self.translate_reserve_peptide, peptide)
					total_reads_to_align += int(future.result()[0])
					to_write_reverse_translation += future.result()[1]

			for peptide, info_peptide in CS_mode.items():
				total_reads_to_align += 1
				sequence = info_peptide[0]
				to_write_reverse_translation+= '>'+peptide+'_'+str(total_reads_to_align)+'\n'+sequence+'\n'

			logging.info('Total Coding Sequences : %d', total_reads_to_align)
			file_to_open = open(path_to_output_folder+'genome_alignments/'+name_exp+'.fastq', 'w')
			file_to_open.write(to_write_reverse_translation)
			file_to_open.close()
			t_2 = time.time()
			total = t_2-t_0
			print ("Total time run function reverse_translation End : %s min" % (total/60.0))
			logging.info('Total time run function reverse_translation to end : %f min', (total/60.0))
		else:
			logging.info('Fasta file with all the coding sequences already exists in the output path : %s --> Skipping this step! ', path_to_output_folder+'genome_alignments/'+name_exp+'.fastq')

	def translate_reserve_peptide(self, peptide):
		sequences = [codon for codon in CODON_TABLE[peptide[0]]]
		count_to_return = 0
		to_print = ''
		for index, amino_acid in enumerate(peptide[1:]):
			to_extend = sequences
			sequences = []
			for index_codon, codon in enumerate(CODON_TABLE[amino_acid]):
				for index_seq,sequence in enumerate(to_extend):
					sequence += codon
					sequences.append(sequence)
					if index == len(peptide[1:])-1 :
						count_to_return += 1
						to_print+= '>'+peptide+'_'+str(count_to_return)+'\n'+sequence+'\n'
		return count_to_return, to_print

