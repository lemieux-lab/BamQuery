import os, logging, time, pickle, multiprocessing, _thread, csv, math, copy, pysam
import pandas as pd
from pathos.multiprocessing import ProcessPool
from collections import Counter
import utils.useful_functions as uf

NUM_WORKERS =  multiprocessing.cpu_count()

path_to_lib = '/'.join(os.path.abspath(__file__).split('/')[:-3])+'/lib/'
annotations_file = path_to_lib+'Info_Transcripts_Annotations.dic'

genome = path_to_lib+'GRCh38.primary_assembly.genome.fa'
genome_index = path_to_lib+'GRCh38.primary_assembly.genome.fa.fai'

__author__ = "Maria Virginia Ruiz Cuevas"


class BiotypeGenomicSearch:

	def __init__(self, peptides_intersected):
		self.peptides_intersected = peptides_intersected


	def get_biotype_from_intersected_transcripts(self):
		
		keys_peptides = []
		keys = []
		transcripts = []
		info_transcripts = []

		with open(annotations_file, 'rb') as fp:
			info_transcripts_dic = pickle.load(fp)

		# Get the information of the transcript from dictionary information
		#{'AAAAPRPAL_chr13:99970439-99970465': {'chr13:99970439-99970465': ['ENST00000267294.4']}}
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
		
		# This gets for the start and the end of the peptide location, where into the transcript the peptide
		# starts or ends, for instance: (5'utr, introns) etc...
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
		
		information_final_biotypes_peptides = self.set_final_transcript_level_biotype(transcripts_information)

		return information_final_biotypes_peptides


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

		presence = ['OutTranscript','OutTranscript']
		for index, point in enumerate([start, end]):
			for key_ in keys:
				present_peptide = self.peptide_in_section(strand, point, info_transcript[key_])
				if present_peptide:
					presence[index] = key_
					break
		
		to_return = [key_peptide, key, transcript, gene_type, transcript_type, presence, info_transcript]
		return to_return

	def peptide_in_section(self, strand, point, info_transcript):

		peptide_into = False
		try:
			r = range(info_transcript[0][0], info_transcript[-1][1]+1)
			if strand == '-':
				r = range(info_transcript[-1][0], info_transcript[0][1]+1)

			if point in r :
				for range_ in info_transcript:
					r = range(range_[0], range_[1]+1)
					peptide_into = point in r
					if peptide_into :
						break
		except IndexError:
			pass
		return peptide_into


	def set_final_transcript_level_biotype(self, transcripts_information):

		information_final_biotypes_peptides = {}

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
							if transcript_level_biotype[0] == 'CDS':
								transcript_level_biotype = self.get_in_frame_out_frame_in_protein(peptide, transcript, info_transcript)
						
						transcript_level_biotype = transcript_level_biotype[0]

						try:
							dic =  information_final_biotypes_peptides[peptide]
							try:
								dic_transcripts = dic[key_peptide]
								dic_transcripts[transcript] = [gene_type, transcript_type, transcript_level_biotype]
							except KeyError:
								dic[key_peptide] = {transcript: [gene_type, transcript_type, transcript_level_biotype]}
						except KeyError:
							dic = {}
							dic[key_peptide] = {transcript: [gene_type, transcript_type, transcript_level_biotype]}
							information_final_biotypes_peptides[peptide] = dic

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
							dic =  information_final_biotypes_peptides[peptide]
							try:
								dic_transcripts = dic[key_peptide]
								dic_transcripts[transcript] = [gene_type, transcript_type, transcript_level_biotype]
							except KeyError:
								dic[key_peptide] = {transcript: [gene_type, transcript_type, transcript_level_biotype]}
						except KeyError:
							dic = {}
							dic[key_peptide] = {transcript: [gene_type, transcript_type, transcript_level_biotype]}
							information_final_biotypes_peptides[peptide] = dic

		return information_final_biotypes_peptides

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

