import os, logging, time, pickle, multiprocessing, _thread, csv, math, copy, pysam
import pandas as pd
from pathos.multiprocessing import ProcessPool
import utils.useful_functions as uf



NUM_WORKERS =  multiprocessing.cpu_count()

path_to_lib = '/'.join(os.path.abspath(__file__).split('/')[:-3])+'/lib/'
annotations_file = path_to_lib+'Info_Transcripts_Annotations.dic'

genome = path_to_lib+'GRCh38.primary_assembly.genome.fa'
genome_index = path_to_lib+'GRCh38.primary_assembly.genome.fa.fai'

__author__ = "Maria Virginia Ruiz Cuevas"

class BiotypeAssignation:

	def __init__(self, path_to_output_folder):
		path_to_lib = '/'.join(os.path.abspath(__file__).split('/')[:-3])+'/lib/'
		self.path_to_output_folder = path_to_output_folder+'res/BED_files/'

		
	def get_information_from_BED_intersection(self):

		file_to_read = self.path_to_output_folder+'intersection_with_annotations.bed'

		self.peptides_intersected = {}

		with open(file_to_read) as f:

			for index, line in enumerate(f):
				splitLine = line.strip().split('\t')
				chr = splitLine[0]
				start = splitLine[1]
				end = splitLine[2]
				key_peptide = splitLine[3]
				strand_peptide = splitLine[4]
				strand_transcript = splitLine[11]
				
				if strand_peptide == strand_transcript :
					split_key = key_peptide.split('_')[1].split('|')
					#chrX	155485562	155485588	TTRPALQEL_chrX:155485562-155485588	-	chrX	HAVANA	transcript	155466540	155494110	.	+	.	gene_id "ENSG00000224533.4"; transcript_id "ENST00000433624.1"; gene_type "antisense"; gene_name "TMLHE-AS1"; transcript_type "antisense"; transcript_name "TMLHE-AS1-002"; level 2; transcript_support_level "1"; tag "basic"; havana_gene "OTTHUMG00000022672.1"; havana_transcript "OTTHUMT00000058815.1";	26
					transcript = ''
					transcript_support_level = ''
					gene_type = ''

					try:
						transcript = line.split(' transcript_id ')[1].split('\"')[1]
					except IndexError:
						pass

					key_peptide_position = chr+':'+start+'-'+end
					
					try:
						self.peptides_intersected[key_peptide][key_peptide_position].append(transcript)
					except KeyError:
						dic = {}
						for key in split_key:
							if 'chr' not in key:
								key = chr+':'+key
							if key == key_peptide_position:
								dic[key] = [transcript]
							else:
								dic[key] = []
						self.peptides_intersected[key_peptide] = dic


	def get_biotype_from_intersected_transcripts(self):
		
		keys_peptides = []
		keys = []
		transcripts = []
		info_transcripts = []

		with open(annotations_file, 'rb') as fp:
			info_transcripts_dic = pickle.load(fp)

		for key_peptide, info_keys in self.peptides_intersected.items():

			for key, intersected_transcripts in info_keys.items():

				for transcript in intersected_transcripts:
					info_transcript = info_transcripts_dic[transcript]
					transcripts.append(transcript)
					info_transcripts.append(info_transcript)

				if len(intersected_transcripts) > 0:
					keys.extend([key]*len(intersected_transcripts))
					keys_peptides.extend([key_peptide]*len(intersected_transcripts))

		self.information_biotypes_peptides = {}
		
		pool = ProcessPool(nodes=NUM_WORKERS)
		results = pool.map(self.biotype_gene_and_transcript_level, keys_peptides, keys, transcripts, info_transcripts)

		transcripts_information = {}
		for res in results:
			key_peptide = res[0]
			key = res[1]  
			transcript = res[2]  
			gene_type = res[3] 
			transcript_type = res[4]  
			presence = res[5]
			info_transcript = res[6] 
			try:
				keys_split_peptides =  self.information_biotypes_peptides[key_peptide]
				try:
					keys_split_peptides[key][transcript] = [gene_type, transcript_type, presence]
				except KeyError:
					keys_split_peptides[key] = {transcript: [gene_type, transcript_type, presence]}
			except KeyError:
				self.information_biotypes_peptides[key_peptide] = {key: {transcript: [gene_type, transcript_type, presence]}}
			
			transcripts_information[transcript] = info_transcript

		#print (self.information_biotypes_peptides)
		self.set_final_transcript_level_biotype(transcripts_information)


	def biotype_gene_and_transcript_level(self, key_peptide, key, transcript, info_transcript):

		gene_type = info_transcript['Info'][5]
		transcript_type = info_transcript['Info'][8]
		peptide = key_peptide.split('_')[0]
		to_return = [key_peptide, key, transcript, gene_type]

		chr = key.split(':')[0]
		start = int(key.split(':')[1].split('-')[0])
		end = int(key.split(':')[1].split('-')[1])
		strand = info_transcript['Info'][3]
		
		if len(info_transcript['CDS']) == 0:
			keys = ['5UTR', '3UTR', 'Introns', 'Exons']
		else:
			keys = ['5UTR', '3UTR', 'Introns', 'CDS']

		presence = ['','']
		for index, point in enumerate([start, end]):
			for key_ in keys:
				present_peptide = self.peptide_in_section(strand, point, info_transcript[key_])
				if present_peptide:
					presence[index] = key_
					break
		
		to_return = [key_peptide, key, transcript, gene_type, transcript_type, presence, info_transcript]
		return to_return

	def peptide_in_section(self, strand, start, info_transcript):

		start_peptide_into = False
		
		try:
			r = range(info_transcript[0][0], info_transcript[-1][1]+1)
			if strand == '-':
				r = range(info_transcript[-1][0], info_transcript[0][1]+1)

			if start in r :
				for range_ in info_transcript:
					r = range(range_[0], range_[1]+1)
					start_peptide_into = start in r
					if start_peptide_into :
						break
		except IndexError:
			pass
		return start_peptide_into


	def set_final_transcript_level_biotype(self, transcripts_information):

		self.information_final_biotypes_peptides = {}

		for key_peptide, positions in self.information_biotypes_peptides.items():
			peptide = key_peptide.split('_')[0]
			position = key_peptide.split('_')[1]

			if '|' not in position:
				for key, transcripts_intersected in positions.items():
					for transcript, biotype in transcripts_intersected.items():
						gene_type = biotype[0]
						transcript_type = biotype[1]
						transcript_level_biotype = list(set(biotype[2]))
						info_transcript = transcripts_information[transcript]

						if len(transcript_level_biotype) > 1 :
							transcript_level_biotype = ['-'.join(biotype[2])]
						else:
							if transcript_level_biotype == 'CDS':
								transcript_level_biotype = self.get_in_frame_out_frame_in_protein(peptide, transcript, info_transcript)
						
						transcript_level_biotype = transcript_level_biotype[0]
						try:
							keys_split_peptides =  self.information_final_biotypes_peptides[key_peptide]
							keys_split_peptides[transcript] = [gene_type, transcript_type, transcript_level_biotype]
						except KeyError:
							self.information_final_biotypes_peptides[key_peptide] = {transcript: [gene_type, transcript_type, transcript_level_biotype]}

			else:
				for key, transcripts_intersected in positions.items():
					for transcript, biotype in transcripts_intersected.items():
						gene_type = biotype[0]
						transcript_type = biotype[1]
						transcript_level_biotype = set(biotype[2])
						info_transcript = transcripts_information[transcript]

						info_transcripts_in_other_positions = []
						for key_2, transcripts_intersected_2 in positions.items():
							if key != key_2:
								try:
									biotype_in_transcript_in_other_position = transcripts_intersected_2[transcript]
									transcript_level_biotype = transcript_level_biotype.union(set(biotype_in_transcript_in_other_position[2]))
									del transcripts_intersected_2[transcript]

								except KeyError:
									pass
						transcript_level_biotype = list(transcript_level_biotype)

						if len(transcript_level_biotype) > 1 :
							transcript_level_biotype = ['-'.join(transcript_level_biotype)]
						else:
							if transcript_level_biotype[0] == 'CDS':
								transcript_level_biotype = self.get_in_frame_out_frame_in_protein(peptide, transcript, info_transcript)

						transcript_level_biotype = transcript_level_biotype[0]	
						try:
							keys_split_peptides =  self.information_final_biotypes_peptides[key_peptide]
							keys_split_peptides[transcript] = [gene_type, transcript_type, transcript_level_biotype]
						except KeyError:
							self.information_final_biotypes_peptides[key_peptide] =  {transcript: [gene_type, transcript_type, transcript_level_biotype]}
		#print (self.information_final_biotypes_peptides)


	def get_in_frame_out_frame_in_protein(self, peptide, transcript, info_transcript):

		chr = info_transcript['Info'][0]
		regions = info_transcript['CDS']
		strand = info_transcript['Info'][3]

		protein = self.get_transcript_and_protein(chr, regions, strand)
				
		if peptide in protein:
			transcript_level = 'In_frame'
		else:
			transcript_level = 'Frameshift'
		return [transcript_level]


	def get_transcript_and_protein(self, chr, regions, strand):

		faFile = pysam.FastaFile(genome, genome_index)

		sequence_transcript = ''
		for cds in regions:
			start_exon = cds[0]
			end_exon = cds[1]
			sequence = faFile.fetch(chr,start_exon-1,end_exon)

			if strand == '-' :
				sequence = uf.reverseComplement(sequence)
				sequence_transcript = sequence_transcript + sequence
				
			elif strand == '+':
				sequence_transcript =  sequence_transcript + sequence

		if chr == 'chrM':
			proteine = uf.translateDNA(sequence_transcript, frame = 'f1', translTable_id='mt')
		else:
			proteine = uf.translateDNA(sequence_transcript, frame = 'f1', translTable_id='default')

		return proteine

		