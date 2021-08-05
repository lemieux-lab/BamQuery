import time, logging, pysam, os
import utils.useful_functions as uf

ATGC = set()
ATGC.add('A')
ATGC.add('T')
ATGC.add('G')
ATGC.add('C')

__author__ = "Maria Virginia Ruiz Cuevas"
__email__ = "maria.virginia.ruiz.cuevas@umontreal.ca"

path_to_lib = '/'.join(os.path.abspath(__file__).split('/')[:-3])+'/lib/'
genomePathFai = path_to_lib + 'GRCh38.primary_assembly.genome.fa.fai'
genomePath = path_to_lib + 'GRCh38.primary_assembly.genome.fa'


class ReadInputFile:

	def __init__(self, path_to_input_folder):
		self.path_to_input_folder = path_to_input_folder

	def treatment_file(self):
		self.peptide_mode = {}
		self.CS_mode = {}
		self.manual_mode = {}
		self.region_quantification_mode = {}
		self.all_mode_peptide = {}
		self.peptides_by_type = {}
		
		peptides_list = self.path_to_input_folder+'peptides.tsv'

		with open(peptides_list) as f:
			for index, line in enumerate(f):
				line = line.strip().split('\t')
				
				if len(line) == 2: 
					peptide = line[0].strip()
					peptide_type = line[1].strip()
					#if not self.evaluate_cs_ntd(peptide_type):
					#	raise Exception("You must provide the peptide type for peptide in the peptide mode instead of the Coding sequence. Otherwise, add the peptide type to evaluate the peptide in CS mode. ", peptide)
					if peptide not in self.all_mode_peptide:
						self.peptide_mode[peptide] = ['','','',peptide_type]
						self.all_mode_peptide[peptide] = peptide_type
											
				elif len(line) == 3:
					if len (line[1].strip()) == 1:
						region = line[0].strip()
						if 'chr' not in region or ':' not in region or '-' not in region:
							raise Exception("Sorry, You must provide an appropriate Genomic Position --> (chr:x-y|z-a) for peptide ", peptide)
						strand = line[1].strip()
						peptide_type = line[2].strip()
						cs, peptide = self.get_local_reference(region, strand)

						if peptide not in self.all_mode_peptide:
							self.region_quantification_mode[peptide] = [cs, region, strand, peptide_type]
							self.all_mode_peptide[peptide] = peptide_type
					else:
						peptide = line[0].strip()
						cs = line[1].strip()
						if self.evaluate_cs_ntd(cs):
							raise Exception("Sorry, You must provide an appropriate Coding Sequence for peptide ", peptide)
						if len(cs)/3 != len(peptide):
							raise Exception("Sorry, Nucleotide sequence doesn\'t correspond to the peptide. ' ", peptide, cs)
						peptide_type = line[2].strip()

						if peptide not in self.all_mode_peptide:
							self.CS_mode[peptide] = [cs,'','',peptide_type]
							self.all_mode_peptide[peptide] = peptide_type

				elif len(line) == 5:
					peptide = line[0].strip()
					cs = line[1].strip()
					if self.evaluate_cs_ntd(cs):
						raise Exception("Sorry, You must provide an appropriate Coding Sequence for peptide ", peptide)
					if len(cs)/3 != len(peptide):
						raise Exception("Sorry, Nucleotide sequence doesn\'t correspond to the peptide. ' ", peptide, cs)
					position = line[2].strip()
					if 'chr' not in position or ':' not in position or '-' not in position:
						raise Exception("Sorry, You must provide an appropriate Genomic Position --> (chr:x-y|z-a) for peptide ", peptide)
					strand = line[3].strip()
					if '+' not in strand and '-' not in strand:
						raise Exception("Sorry, You must provide an appropriate strand, either + (Forward) or - (Backward) for peptide ", peptide)
					peptide_type = line[4].strip()

					if peptide not in self.all_mode_peptide:
						self.manual_mode[peptide] = [cs, position, strand, peptide_type]
						self.all_mode_peptide[peptide] = peptide_type
				else:
					raise Exception("Sorry, You must provide peptides in the appropriate format mode : peptide/CS/manual ")
				
				types = peptide_type.split(';')
				for type_ in types:
					try:
						self.peptides_by_type[type_].append(peptide)
					except KeyError:
						self.peptides_by_type[type_] = [peptide]

		logging.info('Peptides to evaluate in Peptide Mode : %d', len(self.peptide_mode))
		logging.info('Peptides to evaluate in Coding Sequence (CS) Mode : %d', len(self.CS_mode))
		logging.info('Peptides to evaluate in Manual Mode: %d', len(self.manual_mode))
		logging.info('Total Peptides to evaluate : %d', len(self.peptide_mode)+len(self.CS_mode)+len(self.manual_mode))


	def evaluate_cs_ntd(self, cs):
		cs_set = set(cs)
		return len(cs_set - ATGC) != 0

	def get_local_reference(self, region, strand):

		faFile = pysam.FastaFile(genomePath, genomePathFai)
		seqReference = ''
		
		chr = region.split(':')[0]
		regions = region.split(':')[1].split('|')
		for pos in regions:
			start = int(pos.split('-')[0])
			end = int(pos.split('-')[1])
			seqReference += faFile.fetch(chr, start-1,end)

		if strand == '-':
			seqReference = uf.reverseComplement(seqReference)

		if chr=='chrM':
			peptide = uf.translateDNA(seqReference, frame = 'f1', translTable_id='mt')
		else:
			peptide = uf.translateDNA(seqReference, frame = 'f1', translTable_id='default')
		return seqReference, peptide


