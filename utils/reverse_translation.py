import warnings, time, os, gc, pickle
warnings.filterwarnings("ignore")
import billiard as mp

__author__ = "Maria Virginia Ruiz Cuevas"
__email__ = "maria.virginia.ruiz.cuevas@umontreal.ca"

NUM_WORKERS =  int(mp.cpu_count()/2)+2

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
		
		output_message = ''
		self.output_file = path_to_output_folder+'genome_alignments/'+name_exp+'.fastq'
		sam_dic = path_to_output_folder+'genome_alignments/Aligned.out.sam.dic'
		
		exists_fastq = os.path.exists(self.output_file)
		exists_sam_dic = os.path.exists(sam_dic)
		peptide_number = {}
		global q 
		

		if not exists_fastq and not exists_sam_dic:
			t_0 = time.time()
			total_reads_to_align = 0
			to_write_reverse_translation = ''
		
			peptides_list = list(peptide_mode.keys())
			
			#list_of_list_peptides = [peptides_list[x:x+100] for x in range(0, len(peptides_list), 100)]
			
			total_pep = 0
			
			manager = mp.Manager()
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

			q.put('KILL')
			pool.close()
			pool.join()

			self.translate_reserve_peptide_2(peptides_to_re_do)
			
			for peptide, info_peptide in CS_mode.items():
				total_reads_to_align += 1
				sequence = info_peptide[0]
				to_write_reverse_translation+= '>'+peptide+'\n'+sequence+'\n'

			output_message += 'Total Coding Sequences : ' + str(total_reads_to_align)+'\n'
			
			file_to_open = open(self.output_file, 'a')
			file_to_open.write(to_write_reverse_translation)
			file_to_open.close()
			
			with open(path_to_output_folder+'genome_alignments/peptide_number_MCS.dic', 'wb') as handle:
				pickle.dump(peptide_number, handle, protocol=pickle.HIGHEST_PROTOCOL)

			t_2 = time.time()
			total = t_2-t_0

			output_message += 'Total time run function reverse_translation to end : '+ str((total/60.0)) +'min' 
		else:
			path = path_to_output_folder+'genome_alignments/'+name_exp+'.fastq'
			output_message += 'Fasta file with all the coding sequences already exists in the output path : '+path+' --> Skipping this step!'
		return output_message


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
			count_to_return_aux = 0
			
			for index, amino_acid in enumerate(peptide[1:]):
				to_extend = sequences
				sequences = []
				for index_codon, codon in enumerate(CODON_TABLE[amino_acid]):
					for index_seq,sequence in enumerate(to_extend):
						sequence += codon
						sequences.append(sequence)

						if index == len(peptide[1:])-1 :
							count_to_return += 1
							count_to_return_aux += 1
							to_print+= '>'+peptide+'\n'+sequence+'\n'

							if count_to_return_aux == 10000000:
								count_to_return_aux = 0
								sequences = []

		q.put(to_print)
		to_print = ''
		sequences = []
		return count_to_return, peptide

	def translate_reserve_peptide_2(self, peptides):
		
		with open(self.output_file, 'a') as oC :
			for entry in peptides :
				seq = entry.rstrip()
				sequences = [codon for codon in CODON_TABLE[seq[0]]]
				count_to_return = 0
				to_print = ''

				for index, amino_acid in enumerate(seq[1:]):
					to_extend = sequences
					sequences = []

					for codon in CODON_TABLE[amino_acid]:
						for sequence in to_extend:
							sequence += codon
							sequences.append(sequence)

							if index == len(seq[1:])-1 :
								count_to_return += 1
								to_print+= '>'+seq+'\n'+sequence+'\n'
								
								if count_to_return == 1000000:
									oC.write(to_print)
									count_to_return = 0
									to_print = ''
									sequences = []

				if len(to_print) > 0:
					oC.write(to_print)
					to_print = ''


	def listener(self, q):
		
		with open(self.output_file, 'w') as f:
			while 1:
				m = q.get()
				if m == 'KILL':
					break
				f.write(m)
				f.flush()
				m = ''

