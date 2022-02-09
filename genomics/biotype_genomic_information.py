import os, logging, time, pickle, multiprocessing, _thread, csv, math, copy, pysam, time
import pandas as pd
from pathos.multiprocessing import ProcessPool
from collections import Counter
import utils.useful_functions as uf

NUM_WORKERS =  multiprocessing.cpu_count()

path_to_lib = '/'.join(os.path.abspath(__file__).split('/')[:-3])+'/lib/'

__author__ = "Maria Virginia Ruiz Cuevas"


class BiotypeGenomicSearch:

	def __init__(self, peptides_intersected, genome_version, path_to_output_folder):
		self.peptides_intersected = peptides_intersected
		self.path_to_output_folder_alignments = path_to_output_folder+'alignments/'

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

		path = self.path_to_output_folder_alignments +'/alignments_summary_information.pkl'
		self.alignments_summary_information = pd.read_pickle(path)

	def get_biotype_from_intersected_transcripts(self):
		
		t_0 = time.time()

		self.information_biotypes_peptides = {}

		try:
			with open(self.annotations_file, 'rb') as fp:
				info_transcripts_dic = pickle.load(fp)
		except ValueError:
			import pickle5
			with open(self.annotations_file, 'rb') as fp:
				info_transcripts_dic = pickle5.load(fp)

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
						if tsl != '4' or tsl != '5': # Only well supported transcripts are taken into account for biotype calculation. https://m.ensembl.org/info/genome/genebuild/transcript_quality_tags.html
							transcripts_to_search[transcript] = [info_transcript, [key_peptide]]
		
		pool_ = ProcessPool(nodes=NUM_WORKERS)
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
		t_2 = time.time()
		total = t_2-t_0
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
			main_key = key_peptide.split(':')[1]
			chr = main_key.split(':')[0]
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
							#if transcript_type != 'protein_coding' and transcript_type != 'IG_V_gene' and transcript_type != 'TEC' :
							if transcript_type not in ['protein_coding', 'IG_V_gene', 'TEC']:
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

						info_transcripts_in_other_positions = []
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
		
		protein = self.get_transcript_and_protein(chr, regions, strand, len_prot)
		peptide_in_reference = self.alignments_summary_information[(self.alignments_summary_information['Peptide'] == peptide) & (self.alignments_summary_information['Alignment'] == position)]['Peptide in Reference'].values[0] 
		
		if peptide in protein:
			transcript_level = 'In_frame'
		if peptide_in_reference != peptide :
			transcript_level = 'Mutated'
		else:
			transcript_level = 'Frameshift'

		return [transcript_level]


	def get_transcript_and_protein(self, chr, regions, strand, len_prot):

		faFile = pysam.FastaFile(self.genome, self.genome_index)

		sequence_transcript = ''
		check_int = len_prot.is_integer()

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

		if check_int and  '*' not in proteine:
			return proteine
		else:
			if '*' not in proteine:
				return proteine

			if chr == 'chrM':
				proteine = uf.translateDNA(sequence_transcript[1:], frame = 'f1', translTable_id='mt')
			else:
				proteine = uf.translateDNA(sequence_transcript[1:], frame = 'f1', translTable_id='default')
			if '*' in proteine:
				if chr == 'chrM':
					proteine = uf.translateDNA(sequence_transcript[2:], frame = 'f1', translTable_id='mt')
				else:
					proteine = uf.translateDNA(sequence_transcript[2:], frame = 'f1', translTable_id='default')

		return proteine


