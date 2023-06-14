import os, time, pickle, pysam, time
import pandas as pd
from pathos.multiprocessing import ProcessPool
import utils.useful_functions as uf
from Bio import pairwise2

path_to_lib = '/'.join(os.path.abspath(__file__).split('/')[:-3])+'/lib/'

__author__ = "Maria Virginia Ruiz Cuevas"


class BiotypeGenomicSearch:

	def __init__(self, peptides_intersected, genome_version, path_to_output_folder, mouse, threads):
		self.peptides_intersected = peptides_intersected
		self.path_to_output_folder_alignments = path_to_output_folder+'alignments/'
		self.mouse = mouse
		self.threads = threads

		if genome_version == 'v26_88': 
			self.genome_index = path_to_lib + 'genome_versions/genome_v26_88/GRCh38.primary_assembly.genome.fa.fai'
			self.genome = path_to_lib + 'genome_versions/genome_v26_88/GRCh38.primary_assembly.genome.fa'
			self.annotations_file = path_to_lib+'genome_versions/genome_v26_88/Info_Transcripts_Annotations.dic'
		elif genome_version == 'v33_99':
			self.genome_index = path_to_lib + 'genome_versions/genome_v33_99/GRCh38.primary_assembly.genome.fa.fai'
			self.genome = path_to_lib + 'genome_versions/genome_v33_99/GRCh38.primary_assembly.genome.fa'
			self.annotations_file = path_to_lib+'genome_versions/genome_v33_99/Info_Transcripts_Annotations.dic'
		else:
			self.genome_index = path_to_lib + 'genome_versions/genome_v38_104/GRCh38.primary_assembly.genome.fa.fai'
			self.genome = path_to_lib + 'genome_versions/genome_v38_104/GRCh38.primary_assembly.genome.fa'
			self.annotations_file = path_to_lib+'genome_versions/genome_v38_104/Info_Transcripts_Annotations.dic'

		if self.mouse:
			if genome_version == 'M24':
				self.genome_index = path_to_lib + 'genome_versions/genome_mouse_m24/GRCm38.primary_assembly.genome.fa.fai'
				self.genome = path_to_lib + 'genome_versions/genome_mouse_m24/GRCm38.primary_assembly.genome.fa'
				self.annotations_file = path_to_lib+'genome_versions/genome_mouse_m24/Info_Transcripts_Annotations.dic'

			if genome_version == 'M30':
				self.genome_index = path_to_lib + 'genome_versions/genome_mouse_m30/GRCm39.primary_assembly.genome.fa.fai'
				self.genome = path_to_lib + 'genome_versions/genome_mouse_m30/GRCm39.primary_assembly.genome.fa'
				self.annotations_file = path_to_lib+'genome_versions/genome_mouse_m30/Info_Transcripts_Annotations.dic'

		self.translated_prots = {}


	def get_biotype_from_intersected_transcripts(self):
		
		self.information_biotypes_peptides = {}

		with open(self.annotations_file, 'rb') as fp:
			info_transcripts_dic = pickle.load(fp)

		# Get the information of the transcript from dictionary information
		#{'AAAAPRPAL_chr13:99970439-99970465': {'chr13:99970439-99970465': ['ENST00000267294.4']}}
		# NYYPYTITEY	chr17:67552321-67552345|67553610-67553614	AACTATTATCCCTACACAATTACAGAATAC	+	ENST00000581923.5

		transcripts_to_search = {}

		
		for key_peptide, info_keys in self.peptides_intersected.items():

			for key, transcripts_intersected in info_keys.items():
				
				for transcript in transcripts_intersected:
					try:
						transcripts_to_search[transcript][1].append(key_peptide)
					except KeyError:
						info_transcript = info_transcripts_dic[transcript]
						tsl = info_transcript['Info'][9]
						#if tsl != '4' and tsl != '5' : # Only well supported transcripts are taken into account for biotype calculation. https://m.ensembl.org/info/genome/genebuild/transcript_quality_tags.html
						transcripts_to_search[transcript] = [info_transcript, [key_peptide]]
					
		pool_ = ProcessPool(nodes=self.threads)
		results = pool_.map(self.biotype_gene_and_transcript_level, list(transcripts_to_search.values()))
		pool_.close()
		pool_.join()
		pool_.clear()

		# This gets for the start and the end of the peptide location, where into the transcript the peptide
		# starts or ends, for instance: (5'utr, introns) etc...
		
		transcripts_information = {}
		for result in results:
			for res in result:
				key_peptide = res[0]
				key = res[1]  
				transcript = res[2]  
				gene_type = res[3] 
				transcript_type = res[4]  
				presence = res[5]
				info_transcript = transcripts_to_search[transcript][0]

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
		self.translated_prots = {}
		return information_final_biotypes_peptides


	def biotype_gene_and_transcript_level(self, info_transcript_intersected ):

		info_transcript = info_transcript_intersected[0]
		key_peptides = info_transcript_intersected[1]

		gene_type = info_transcript['Info'][5]
		transcript_type = info_transcript['Info'][8]
		transcript = info_transcript['Info'][7]
		to_return = []

		
		for key_peptide in key_peptides:
			peptide = key_peptide.split('_')[0]
			position = key_peptide.split('_')[1]
			main_key = key_peptide.split(':')[1].split('_')[0]
			keys = main_key.split('|')
			
			for index, key in enumerate(keys):
				
				start = int(key.split('-')[0])
				end = int(key.split('-')[1])
				strand = info_transcript['Info'][3]

				
				if len(info_transcript['CDS']) == 0:
					keys_bio = ['Introns', 'Exons']
				else:
					keys_bio = ['5UTR', '3UTR', 'Introns', 'CDS']

				presence = ['Intergenic','Intergenic']
				
				for index, point in enumerate([start, end]):
					
					for key_ in keys_bio:

						present_peptide = self.peptide_in_section(strand, point, info_transcript[key_])
						
						if present_peptide:
							if transcript_type not in ['protein_coding', 'IG_C_gene', 'IG_D_gene', 'IG_J_gene', 'IG_V_gene', 'TR_C_gene', 'TR_D_gene', 'TR_J_gene', 'TR_V_gene']:
								if key_ == 'Exons' or key_ == '5UTR' or key_ == '3UTR':
									key_ = 'Non_coding Exons'

							presence[index] = key_
							break
				
				to_return.append([key_peptide, key, transcript, gene_type, transcript_type, presence])

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
			mcs = key_peptide.split('_')[2]
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
								transcript_level_biotype = self.get_in_frame_out_frame_in_protein(peptide, transcript, info_transcript, position)

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
						transcript_level_biotype = list(set(biotype[2]))
						info_transcript = transcripts_information[transcript]

						for key_2, transcripts_intersected_2 in positions.items():
							if key != key_2:
								try:
									biotype_in_transcript_in_other_position = transcripts_intersected_2[transcript]
									transcript_level_biotype.extend(list(set(biotype_in_transcript_in_other_position[2])))
									del transcripts_intersected_2[transcript]

								except KeyError:
									pass

						transcript_level_biotype_set = set(transcript_level_biotype)
						
						if len(transcript_level_biotype_set) > 1 :
							transcript_level_biotype = ['-'.join(transcript_level_biotype)]
						else:
							if transcript_level_biotype[0] == 'CDS':
								transcript_level_biotype = self.get_in_frame_out_frame_in_protein(peptide, transcript, info_transcript, position)
								
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


	def get_in_frame_out_frame_in_protein(self, peptide, transcript, info_transcript, position):

		chr = info_transcript['Info'][0]
		regions = info_transcript['CDS']
		strand = info_transcript['Info'][3]
		len_prot = info_transcript['Info'][13]
		
		try:
			proteins = self.translated_prots[transcript]
		except KeyError:
			proteins = self.get_transcript_and_protein(chr, regions, strand)
			self.translated_prots[transcript] = proteins
		
		alignments = pairwise2.align.localms(proteins[0], peptide, 1, -1, -5, -.1)
		alignment_score = alignments[0].score
		percentage_similarity = (alignment_score / len(peptide)) 
		print (peptide, percentage_similarity, proteins[0])

		if percentage_similarity >= 0.5:
			transcript_level = 'In_frame'
		else:
			alignments_1 = pairwise2.align.localms(proteins[1], peptide, 1, -1, -5, -.1)
			alignment_score_1 = alignments_1[0].score
			percentage_similarity_1 = (alignment_score_1 / len(peptide)) 
			alignments_2 = pairwise2.align.localms(proteins[2], peptide, 1, -1, -5, -.1)
			alignment_score_2 = alignments_2[0].score
			percentage_similarity_2 = (alignment_score_2 / len(peptide)) 
			
			if percentage_similarity_1 >= 0.5 or percentage_similarity_2 >= 0.5 : 
				transcript_level = 'Frameshift'
			else:
				transcript_level = 'CDS'
	
		return [transcript_level]


	def get_transcript_and_protein(self, chr, regions, strand):

		faFile = pysam.FastaFile(self.genome, self.genome_index)
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
			protein_1 = uf.translateDNA(sequence_transcript, frame = 'f1', translTable_id='mt')
		else:
			protein_1 = uf.translateDNA(sequence_transcript, frame = 'f1', translTable_id='default')
		
		if chr == 'chrM':
			protein_2 = uf.translateDNA(sequence_transcript, frame = 'f2', translTable_id='mt')
		else:
			protein_2 = uf.translateDNA(sequence_transcript, frame = 'f2', translTable_id='default')

		if chr == 'chrM':
			protein_3 = uf.translateDNA(sequence_transcript, frame = 'f3', translTable_id='mt')
		else:
			protein_3 = uf.translateDNA(sequence_transcript, frame = 'f3', translTable_id='default')

		return [protein_1, protein_2, protein_3]


