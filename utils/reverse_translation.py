import warnings, time, logging, os, gc, pickle
warnings.filterwarnings("ignore")
import concurrent.futures
from pathos.multiprocessing import ProcessPool
import multiprocessing as mp

__author__ = "Maria Virginia Ruiz Cuevas"
__email__ = "maria.virginia.ruiz.cuevas@umontreal.ca"

NUM_WORKERS =  mp.cpu_count()

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


class ReverseTranslation:

	def __init__(self):
		pass

	def reverse_translation(self, peptide_mode, CS_mode, path_to_output_folder, name_exp):
		
		self.output_file = path_to_output_folder+'genome_alignments/'+name_exp+'.fastq'
		exists = os.path.exists(self.output_file)
		exists_light = os.path.exists(path_to_output_folder+'alignments/Alignments_information_light.dic')
		peptide_number = {}


		if not exists and not exists_light:
			t_0 = time.time()
			total_reads_to_align = 0
			to_write_reverse_translation = ''
		
			peptides_list = list(peptide_mode.keys())
			list_of_list_peptides = [peptides_list[x:x+100] for x in range(0, len(peptides_list), 100)]
			
			total_pep = 0
			
			manager = mp.Manager()
			global q   
			q = manager.Queue() 
			pool = mp.Pool(NUM_WORKERS)
			watcher = pool.apply_async(self.listener, (q,))
			
			results = pool.map(self.translate_reserve_peptide, peptides_list) 
			
			peptides_to_re_do = []

			for res in results:
				count_to_return = int(res[0])
				peptide = res[1]

				total_pep += 1
				total_reads_to_align += count_to_return
				peptide_number[peptide] = count_to_return

				if count_to_return >= 30000000:
					peptides_to_re_do.append(peptide)

			q.put('kill')
			pool.close()
			pool.join()

			self.translate_reserve_peptide_2(peptides_to_re_do)
			
			for peptide, info_peptide in CS_mode.items():
				total_reads_to_align += 1
				sequence = info_peptide[0]
				to_write_reverse_translation+= '>'+peptide+'\n'+sequence+'\n'

			logging.info('Total Coding Sequences : %d', total_reads_to_align)
			
			file_to_open = open(self.output_file, 'a')
			file_to_open.write(to_write_reverse_translation)
			file_to_open.close()
			
			with open(path_to_output_folder+'genome_alignments/peptide_number_MCS.dic', 'wb') as handle:
				pickle.dump(peptide_number, handle, protocol=pickle.HIGHEST_PROTOCOL)

			t_2 = time.time()
			total = t_2-t_0

			logging.info('Total time run function reverse_translation to end : %f min', (total/60.0))
		else:
			if exists_light:
				logging.info('Skipping : Reverse Translation! ')
			else:
				logging.info('Fasta file with all the coding sequences already exists in the output path : %s --> Skipping this step! ', path_to_output_folder+'genome_alignments/'+name_exp+'.fastq')


	def translate_reserve_peptide(self, peptide):
		list_aa = list(peptide)
		count_to_return = 1
		to_print = ''

		for aa in list_aa:
			count_to_return *= len(CODON_TABLE[aa])
		if count_to_return >= 30000000:
			return count_to_return, peptide

		else:
			sequences = [codon for codon in CODON_TABLE[peptide[0]]]
			count_to_return = 0
			
			for index, amino_acid in enumerate(peptide[1:]):
				to_extend = sequences
				sequences = []
				for index_codon, codon in enumerate(CODON_TABLE[amino_acid]):
					for index_seq,sequence in enumerate(to_extend):
						sequence += codon
						sequences.append(sequence)
						if index == len(peptide[1:])-1 :
							count_to_return += 1
							to_print+= '>'+peptide+'\n'+sequence+'\n'

			q.put(to_print)
			to_print = ''

		return count_to_return, peptide

	def translate_reserve_peptide_2(self, peptides):
		file_to_open = open(self.output_file, 'a')

		for peptide in peptides:
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
							to_print += '>'+peptide+'\n'+sequence+'\n'
							if count_to_return == 100000000:
								file_to_open.write(to_print)
								file_to_open.flush()
								to_print = ''
								count_to_return = 0
			file_to_open.write(to_print)
			file_to_open.flush()


		file_to_open.close()


	def listener(self, q):
		
		with open(self.output_file, 'w') as f:
			while 1:
				m = q.get()
				if m == 'kill':
					break
				f.write(m)
				f.flush()
				m = ''

