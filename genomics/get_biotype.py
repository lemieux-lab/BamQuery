import os, logging, time, pickle, multiprocessing, _thread, csv, math, copy, pysam
import pandas as pd
from pathos.multiprocessing import ProcessPool
import utils.useful_functions as uf
import plotting.plots as plots
from collections import Counter


NUM_WORKERS =  multiprocessing.cpu_count()

path_to_lib = '/'.join(os.path.abspath(__file__).split('/')[:-3])+'/lib/'
annotations_file = path_to_lib+'Info_Transcripts_Annotations.dic'

genome = path_to_lib+'GRCh38.primary_assembly.genome.fa'
genome_index = path_to_lib+'GRCh38.primary_assembly.genome.fa.fai'

__author__ = "Maria Virginia Ruiz Cuevas"

class BiotypeAssignation:

	def __init__(self, path_to_output_folder, name_exp, mode):
		path_to_lib = '/'.join(os.path.abspath(__file__).split('/')[:-3])+'/lib/'
		self.path_to_output_folder = path_to_output_folder
		self.path_to_output_folder_bed_files = path_to_output_folder+'res/BED_files/'
		self.name_exp = name_exp
		self.mode = mode
		with open(path_to_lib+'ERE_info.dic', 'rb') as handle:
			self.ere_info = pickle.load(handle)
		

	def get_information_from_BED_intersection(self):

		file_to_read = self.path_to_output_folder_bed_files+'intersection_with_annotated_transcripts.bed'

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
					# chr21	42127754	42127780	FINDSIVY_chr21:42127754-42127780	+	chr21	HAVANA	exon	42127678	42127831	.	+	.	gene_id "ENSG00000177398.18"; transcript_id "ENST00000484712.1"; gene_type "protein_coding"; gene_name "UMODL1"; transcript_type "processed_transcript"; transcript_name "UMODL1-013"; exon_number 1; exon_id "ENSE00001920427.1"; level 2; transcript_support_level "1"; havana_gene "OTTHUMG00000086788.3"; havana_transcript "OTTHUMT00000195303.1";	26
					transcript = ''
					transcript_support_level = ''
					gene_type = ''
					peptide = key_peptide.split('_')[0]
					intersection_number = int(splitLine[-1])

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


	def get_information_from_BED_intersection_ERE(self):

		file_to_read = self.path_to_output_folder_bed_files+'intersection_with_annotated_EREs.bed'

		self.peptides_intersected_ere = {}
		# chrY	17666734	17666760	TTRPALQEL_chrY:17666734-17666760	-	chrY	17664318	17668712	L1P2	14950	-	26
		with open(file_to_read) as f:

			for index, line in enumerate(f):
				splitLine = line.strip().split('\t')
				chr = splitLine[0]
				start = splitLine[1]
				end = splitLine[2]
				key_peptide = splitLine[3]
				strand_peptide = splitLine[4]
				strand_transcript = splitLine[10]
				
				if strand_peptide == strand_transcript :
					split_key = key_peptide.split('_')[1].split('|')
					#chrY	17664318	17668712	L1P2	14950	-	26
					repName = splitLine[8]
					intersection_number = int(splitLine[11])
					peptide = key_peptide.split('_')[0]
					key_peptide_position = key_peptide.split('_')[1]+'_'+strand_peptide
					
					try:
						peptide_info = self.peptides_intersected_ere[peptide]
						try:
							repName_set = peptide_info[key_peptide_position]
							repName_set.add(repName)
						except KeyError:
							repName_set = set()
							repName_set.add(repName)
							peptide_info[key_peptide_position] = repName_set
					except KeyError:
						repName_set = set()
						repName_set.add(repName)
						self.peptides_intersected_ere[peptide] = {key_peptide_position : repName_set}

	

	def get_biotype_from_EREs(self, info_peptide_alignments, peptides_by_type):
		
		data = []
		data_final = []

		biotypes_peptides_rna_ERE = {}
		biotypes_peptides_ribo_ERE = {}

		biotypes_peptides_rna_ere_type_peptide = {}
		biotypes_peptides_ribo_ere_type_peptide = {}

		biotypes_peptides_ere_type_peptide = {}
		biotypes_peptides_ere_all = {}

		for type_peptide, peptides in peptides_by_type.items():
			
			biotypes_peptides_rna_ERE[type_peptide] = {}
			biotypes_peptides_ribo_ERE[type_peptide] = {}
			ere_biotype_type_peptide_rna = []
			ere_biotype_type_peptide_ribo = []

			biotypes_peptides_ere_type_peptide[type_peptide] = {}

			for peptide in peptides:
				
				ere_biotype_peptide = []
				ere_biotype_family_peptide = []
				total_count_rna = 0
				total_count_ribo = 0

				try:
					info_peptide = self.peptides_intersected_ere[peptide]
					
					for key_peptide in info_peptide:
						
						repName = list(info_peptide[key_peptide])
						alignment = key_peptide.split('_')[0]
						strand = key_peptide.split('_')[1]

						alignments = info_peptide_alignments[peptide]
						count_rna = 0
						count_ribo = 0
						for alig in alignments:
							if alignment == alig[0]:
								count_rna = alig[1]
								count_ribo = alig[2]
								break

						total_count_ribo += count_ribo
						total_count_rna += count_rna

						for name in repName:
							repClass = self.ere_info[name][0]
							repFamily = self.ere_info[name][1]

							ere_biotype_peptide.append(repClass)
							ere_biotype_family_peptide.append(repFamily)

							if count_rna > 0:
								ere_biotype_type_peptide_rna.append(repClass)

							if count_ribo > 0:
								ere_biotype_type_peptide_ribo.append(repClass)

							to_add = [type_peptide, peptide, alignment, strand, name, repClass, repFamily, count_rna, count_ribo]
							data.append(to_add)
				except KeyError:
					pass

				count_ere_biotypes = Counter(ere_biotype_peptide)
				total_ere_biotype = len(ere_biotype_peptide)
				type_ere_peptide_biotype = ''

				for ere_biotype, v in count_ere_biotypes.items():
					value =  v / (total_ere_biotype*1.0)

					try:
						biotypes_peptides_ere_type_peptide[type_peptide][ere_biotype] += value
					except KeyError:
						biotypes_peptides_ere_type_peptide[type_peptide][ere_biotype] = value
					
					try:
						biotypes_peptides_ere_all[ere_biotype] += value
					except KeyError:
						biotypes_peptides_ere_all[ere_biotype] = value

					type_ere_peptide_biotype += ere_biotype+': '+str(round(value*100,2))+'% - ' 
				
				if type_ere_peptide_biotype == '':
					type_ere_peptide_biotype = 'NA -'

				count_ere_biotypes = Counter(ere_biotype_family_peptide)
				total_ere_biotype = len(ere_biotype_family_peptide)
				type_ere_biotype_family_peptide = ''

				for ere_biotype, v in count_ere_biotypes.items():
					value =  v / (total_ere_biotype*1.0)
					type_ere_biotype_family_peptide += ere_biotype+': '+str(round(value*100,2))+'% - ' 
				
				if type_ere_biotype_family_peptide == '':
					type_ere_biotype_family_peptide = 'NA -'

				to_add = [type_peptide, peptide, type_ere_peptide_biotype[:-2], type_ere_biotype_family_peptide[:-2], total_count_rna, total_count_ribo]
				data_final.append(to_add)


			labels = list(biotypes_peptides_ere_type_peptide[type_peptide].keys())
			sizes = list(biotypes_peptides_ere_type_peptide[type_peptide].values())
			title = 'ERE biotype '+type_peptide
			if sum(sizes) > 0:
				plots.plot_pie_ere(title, labels, sizes, self.path_to_output_folder+'plots/biotypes/ERE_annotation/global_biotypes/', self.name_exp, type_peptide)

			count_ere_biotypes = Counter(ere_biotype_type_peptide_rna)
			total_ere_biotype = len(ere_biotype_type_peptide_rna)
			
			for ere_biotype, v in count_ere_biotypes.items():
				value =  v / (total_ere_biotype*1.0)

				try:
					biotypes_peptides_rna_ERE[type_peptide][ere_biotype] += value
				except KeyError:
					biotypes_peptides_rna_ERE[type_peptide][ere_biotype] = value

				try:
					biotypes_peptides_rna_ere_type_peptide[ere_biotype] += value
				except KeyError:
					biotypes_peptides_rna_ere_type_peptide[ere_biotype] = value

			labels = list(biotypes_peptides_rna_ERE[type_peptide].keys())
			sizes = list(biotypes_peptides_rna_ERE[type_peptide].values())
			title = 'ERE biotype '+type_peptide +' rna '
			if sum(sizes) > 0:
				plots.plot_pie_ere(title, labels, sizes, self.path_to_output_folder+'plots/biotypes/ERE_annotation/with_counts/', self.name_exp, 'rna_'+type_peptide)

			count_ere_biotypes = Counter(ere_biotype_type_peptide_ribo)
			total_ere_biotype = len(ere_biotype_type_peptide_ribo)
			
			for ere_biotype, v in count_ere_biotypes.items():
				value =  v / (total_ere_biotype*1.0)

				try:
					biotypes_peptides_ribo_ERE[type_peptide][ere_biotype] += value
				except KeyError:
					biotypes_peptides_ribo_ERE[type_peptide][ere_biotype] = value

				try:
					biotypes_peptides_ribo_ere_type_peptide[ere_biotype] += value
				except KeyError:
					biotypes_peptides_ribo_ere_type_peptide[ere_biotype] = value
			
			labels = list(biotypes_peptides_ribo_ERE[type_peptide].keys())
			sizes = list(biotypes_peptides_ribo_ERE[type_peptide].values())
			title = 'ERE biotype '+type_peptide +' ribo '
			if sum(sizes) > 0:
				plots.plot_pie_ere(title, labels, sizes, self.path_to_output_folder+'plots/biotypes/ERE_annotation/with_counts/', self.name_exp, 'ribo_'+type_peptide)

		labels = list(biotypes_peptides_ribo_ere_type_peptide.keys())
		sizes = list(biotypes_peptides_ribo_ere_type_peptide.values())
		title = 'ERE biotype All peptides ribo '
		if sum(sizes) > 0:
			plots.plot_pie_ere(title, labels, sizes, self.path_to_output_folder+'plots/biotypes/ERE_annotation/with_counts/', self.name_exp, 'ribo_All_peptides')

		labels = list(biotypes_peptides_rna_ere_type_peptide.keys())
		sizes = list(biotypes_peptides_rna_ere_type_peptide.values())
		title = 'ERE biotype All peptides ribo '
		if sum(sizes) > 0:
			plots.plot_pie_ere(title, labels, sizes, self.path_to_output_folder+'plots/biotypes/ERE_annotation/with_counts/', self.name_exp, 'rna_All_peptides')


		labels = list(biotypes_peptides_ere_all.keys())
		sizes = list(biotypes_peptides_ere_all.values())
		title = 'ERE biotype All peptides'
		if sum(sizes) > 0:
			plots.plot_pie_ere(title, labels, sizes, self.path_to_output_folder+'plots/biotypes/ERE_annotation/global_biotypes/', self.name_exp, 'All_peptides')


		self.df4 = pd.DataFrame(data, columns=['Peptide Type', 'Peptide','Alignment', 'Strand', 'ERE name', 'ERE class', 'ERE family', 'reads count RNA', 'reads count Ribo'])
		self.df5 = pd.DataFrame(data_final, columns=['Peptide Type', 'Peptide','ERE class','ERE family', 'Total reads count RNA', 'Total reads count Ribo'])


	def get_biotype_from_intersected_transcripts(self):
		
		keys_peptides = []
		keys = []
		transcripts = []
		info_transcripts = []

		with open(annotations_file, 'rb') as fp:
			info_transcripts_dic = pickle.load(fp)

		# Get the information of the transcript from dictionary information
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
							if transcript_level_biotype[0] == 'CDS':
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


	def prepare_info_to_draw_biotypes(self, info_peptide_alignments, peptides_by_type):

		biotypes_peptides_rna = {}
		biotypes_peptides_ribo = {}
		biotypes_peptides_rna_all_peptides = {'gene_level_biotype':{}, 'transcript_level_biotype':{}, 'genomic_position_biotype':{}}
		biotypes_peptides_ribo_all_peptides = {'gene_level_biotype':{}, 'transcript_level_biotype':{}, 'genomic_position_biotype':{}}
		biotypes_peptides_all_peptides = {'gene_level_biotype':{}, 'transcript_level_biotype':{}, 'genomic_position_biotype':{}}

		biotypes_peptides_global = {'gene_level_biotype':{}, 'transcript_level_biotype':{}, 'genomic_position_biotype':{}}
		
		data = []
		data_resume = []
		data_final = []
		
		for type_peptide, peptides in peptides_by_type.items():
			biotypes_peptides_rna[type_peptide] = {'gene_level_biotype':{}, 'transcript_level_biotype':{}, 'genomic_position_biotype':{}}
			biotypes_peptides_ribo[type_peptide] = {'gene_level_biotype':{}, 'transcript_level_biotype':{}, 'genomic_position_biotype':{}}

			biotypes_peptides_global[type_peptide] = {'gene_level_biotype':{}, 'transcript_level_biotype':{}, 'genomic_position_biotype':{}}

			gene_level_biotypes_type_peptide_rna = []
			transcript_level_biotypes_type_peptide_rna = []
			genomic_position_biotypes_type_peptide_rna = []

			gene_level_biotypes_type_peptide_ribo = []
			transcript_level_biotypes_type_peptide_ribo = []
			genomic_position_biotypes_type_peptide_ribo = []

			for peptide in peptides:

				total_count_rna = 0
				total_count_ribo = 0
				gene_level_biotypes_peptide = []
				transcript_level_biotypes_peptide = []
				genomic_position_biotypes_peptide = []

				try:
					alignments = info_peptide_alignments[peptide]
					count_rna = 0
					count_ribo = 0

					for alignment in alignments:
						key = peptide+'_'+alignment[0]
						count_rna = alignment[1]
						count_ribo = alignment[2]
						strand = alignment[3]
						alignment = alignment[0]

						total_count_rna += count_rna
						total_count_ribo += count_ribo
						
						try:
							transcripts_intersected = self.information_final_biotypes_peptides[key]

							gene_level_biotypes = []
							transcript_level_biotypes = []
							genomic_position_biotypes = []

							for position, transcript in transcripts_intersected.items():
								gene_level_biotype = transcript[0]
								transcript_level_biotype = transcript[1]
								genomic_position_biotype = transcript[2]
								
								gene_level_biotypes.append(gene_level_biotype)
								transcript_level_biotypes.append(transcript_level_biotype)
								genomic_position_biotypes.append(genomic_position_biotype)

								to_add = [type_peptide, peptide, alignment, strand, position, gene_level_biotype, transcript_level_biotype, genomic_position_biotype, count_rna, count_ribo]
								
								data.append(to_add)

							count_gene_level_biotypes = Counter(gene_level_biotypes)
							count_transcript_level_biotypes = Counter(transcript_level_biotypes)
							count_genomic_position_biotypes = Counter(genomic_position_biotypes)

							
							total_gene_level_biotype = len(gene_level_biotypes)
							type_gene_level_biotype = ''
							
							for gene_level_biotype, v in count_gene_level_biotypes.items():
								gene_level_biotypes_peptide.append(gene_level_biotype)
								
								value =  v / (total_gene_level_biotype*1.0)

								if count_rna > 0:
									gene_level_biotypes_type_peptide_rna.append(gene_level_biotype)
									try:
										biotypes_peptides_rna[type_peptide]['gene_level_biotype'][gene_level_biotype] += value
									except KeyError:
										biotypes_peptides_rna[type_peptide]['gene_level_biotype'][gene_level_biotype] = value
								
								if count_ribo > 0:
									gene_level_biotypes_type_peptide_ribo.append(gene_level_biotype)
									try:
										biotypes_peptides_ribo[type_peptide]['gene_level_biotype'][gene_level_biotype] += value
									except KeyError:
										biotypes_peptides_ribo[type_peptide]['gene_level_biotype'][gene_level_biotype] = value

								try:
									biotypes_peptides_global[type_peptide]['gene_level_biotype'][gene_level_biotype] += value
								except KeyError:
									biotypes_peptides_global[type_peptide]['gene_level_biotype'][gene_level_biotype] = value

								type_gene_level_biotype += gene_level_biotype+': '+str(round(value*100,2))+'% - ' 
							

							total_transcript_level_biotype = len(transcript_level_biotypes)
							type_transcript_level_biotype = ''
							
							for transcript_level_biotype, v in count_transcript_level_biotypes.items():
								transcript_level_biotypes_peptide.append(transcript_level_biotype)
								
								value =  v / (total_transcript_level_biotype*1.0)
								if count_rna > 0:
									transcript_level_biotypes_type_peptide_rna.append(transcript_level_biotype)
									try:
										biotypes_peptides_rna[type_peptide]['transcript_level_biotype'][transcript_level_biotype] += value
									except KeyError:
										biotypes_peptides_rna[type_peptide]['transcript_level_biotype'][transcript_level_biotype] = value
								if count_ribo > 0:
									transcript_level_biotypes_type_peptide_ribo.append(transcript_level_biotype)
									try:
										biotypes_peptides_ribo[type_peptide]['transcript_level_biotype'][transcript_level_biotype] += value
									except KeyError:
										biotypes_peptides_ribo[type_peptide]['transcript_level_biotype'][transcript_level_biotype] = value

								try:
									biotypes_peptides_global[type_peptide]['transcript_level_biotype'][transcript_level_biotype] += value
								except KeyError:
									biotypes_peptides_global[type_peptide]['transcript_level_biotype'][transcript_level_biotype] = value

								type_transcript_level_biotype += transcript_level_biotype+': '+str(round(value*100,2))+'% - ' 
							
							total_genomic_position_biotype = len(genomic_position_biotypes)
							type_genomic_position_biotype = ''
							
							for genomic_position_biotype, v in count_genomic_position_biotypes.items():
								genomic_position_biotypes_peptide.append(genomic_position_biotype)
								
								value =  v / (total_genomic_position_biotype*1.0)
								if count_rna > 0:
									genomic_position_biotypes_type_peptide_rna.append(genomic_position_biotype)
									try:
										biotypes_peptides_rna[type_peptide]['genomic_position_biotype'][genomic_position_biotype] += value
									except KeyError:
										biotypes_peptides_rna[type_peptide]['genomic_position_biotype'][genomic_position_biotype] = value
								if count_ribo > 0:
									genomic_position_biotypes_type_peptide_ribo.append(genomic_position_biotype)
									try:
										biotypes_peptides_ribo[type_peptide]['genomic_position_biotype'][genomic_position_biotype] += value
									except KeyError:
										biotypes_peptides_ribo[type_peptide]['genomic_position_biotype'][genomic_position_biotype] = value

								try:
									biotypes_peptides_global[type_peptide]['genomic_position_biotype'][genomic_position_biotype] += value
								except KeyError:
									biotypes_peptides_global[type_peptide]['genomic_position_biotype'][genomic_position_biotype] = value


								type_genomic_position_biotype += genomic_position_biotype+': '+str(round(value*100,2))+'% - ' 
							
							to_add = [type_peptide, peptide, alignment, strand, type_gene_level_biotype[:-2], type_transcript_level_biotype[:-2], type_genomic_position_biotype[:-2], count_rna, count_ribo]
							
							data_resume.append(to_add)
	
						except KeyError:

							biotype_type = 'Intergenic'
							gene_level_biotype = 'Intergenic'
							transcript_level_biotype = 'Intergenic'
							genomic_position_biotype = biotype_type

							gene_level_biotypes_peptide.append(gene_level_biotype)
							transcript_level_biotypes_peptide.append(transcript_level_biotype)
							genomic_position_biotypes_peptide.append(genomic_position_biotype)

							to_add = [type_peptide, peptide, alignment, strand, 'No Annotation', gene_level_biotype, transcript_level_biotype, genomic_position_biotype, count_rna, count_ribo]
							data.append(to_add)
							
							to_add = [type_peptide, peptide, alignment, strand, gene_level_biotype, transcript_level_biotype, genomic_position_biotype, count_rna, count_ribo]
							data_resume.append(to_add)

							if count_rna > 0:
								gene_level_biotypes_type_peptide_rna.append(gene_level_biotype)
								transcript_level_biotypes_type_peptide_rna.append(transcript_level_biotype)
								genomic_position_biotypes_type_peptide_rna.append(genomic_position_biotype)

								try:
									biotypes_peptides_rna[type_peptide]['gene_level_biotype'][gene_level_biotype] += 1
								except KeyError:
									biotypes_peptides_rna[type_peptide]['gene_level_biotype'][gene_level_biotype] = 1

								try:
									biotypes_peptides_rna[type_peptide]['transcript_level_biotype'][transcript_level_biotype] += 1
								except KeyError:
									biotypes_peptides_rna[type_peptide]['transcript_level_biotype'][transcript_level_biotype] = 1

								try:
									biotypes_peptides_rna[type_peptide]['genomic_position_biotype'][genomic_position_biotype] += 1
								except KeyError:
									biotypes_peptides_rna[type_peptide]['genomic_position_biotype'][genomic_position_biotype] = 1

							try:
								biotypes_peptides_global[type_peptide]['gene_level_biotype'][gene_level_biotype] += 1
							except KeyError:
								biotypes_peptides_global[type_peptide]['gene_level_biotype'][gene_level_biotype] = 1

							try:
								biotypes_peptides_global[type_peptide]['transcript_level_biotype'][transcript_level_biotype] += 1
							except KeyError:
								biotypes_peptides_global[type_peptide]['transcript_level_biotype'][transcript_level_biotype] = 1

							try:
								biotypes_peptides_global[type_peptide]['genomic_position_biotype'][genomic_position_biotype] += 1
							except KeyError:
								biotypes_peptides_global[type_peptide]['genomic_position_biotype'][genomic_position_biotype] = 1

							if count_ribo > 0:
								gene_level_biotypes_type_peptide_ribo.append(gene_level_biotype)
								transcript_level_biotypes_type_peptide_ribo.append(transcript_level_biotype)
								genomic_position_biotypes_type_peptide_ribo.append(genomic_position_biotype)

								try:
									biotypes_peptides_ribo[type_peptide]['gene_level_biotype'][gene_level_biotype] += 1
								except KeyError:
									biotypes_peptides_ribo[type_peptide]['gene_level_biotype'][gene_level_biotype] = 1

								try:
									biotypes_peptides_ribo[type_peptide]['transcript_level_biotype'][transcript_level_biotype] += 1
								except KeyError:
									biotypes_peptides_ribo[type_peptide]['transcript_level_biotype'][transcript_level_biotype] = 1

								try:
									biotypes_peptides_ribo[type_peptide]['genomic_position_biotype'][genomic_position_biotype] += 1
								except KeyError:
									biotypes_peptides_ribo[type_peptide]['genomic_position_biotype'][genomic_position_biotype] = 1
	
				except KeyError:
					pass

				# Peptides
				count_gene_level_biotypes = Counter(gene_level_biotypes_peptide)
				count_transcript_level_biotypes = Counter(transcript_level_biotypes_peptide)
				count_genomic_position_biotypes = Counter(genomic_position_biotypes_peptide)
				
				total_gene_level_biotype = len(gene_level_biotypes_peptide)
				type_gene_level_biotype = ''
				
				for gene_level_biotype, v in count_gene_level_biotypes.items():
					value =  v / (total_gene_level_biotype*1.0)
					try:
						biotypes_peptides_all_peptides['gene_level_biotype'][gene_level_biotype] += value
					except KeyError:
						biotypes_peptides_all_peptides['gene_level_biotype'][gene_level_biotype] = value

					type_gene_level_biotype += gene_level_biotype+': '+str(round(value*100,2))+'% - ' 
					
				total_transcript_level_biotype = len(transcript_level_biotypes_peptide)
				type_transcript_level_biotype = ''
				
				for transcript_level_biotype, v in count_transcript_level_biotypes.items():
					value =  v / (total_transcript_level_biotype*1.0)
					try:
						biotypes_peptides_all_peptides['transcript_level_biotype'][transcript_level_biotype] += value
					except KeyError:
						biotypes_peptides_all_peptides['transcript_level_biotype'][transcript_level_biotype] = value

					type_transcript_level_biotype += transcript_level_biotype+': '+str(round(value*100,2))+'% - ' 
					
				total_genomic_position_biotype = len(genomic_position_biotypes_peptide)
				type_genomic_position_biotype = ''
				
				for genomic_position_biotype, v in count_genomic_position_biotypes.items():
					value =  v / (total_genomic_position_biotype*1.0)
					try:
						biotypes_peptides_all_peptides['genomic_position_biotype'][genomic_position_biotype] += value
					except KeyError:
						biotypes_peptides_all_peptides['genomic_position_biotype'][genomic_position_biotype] = value

					type_genomic_position_biotype += genomic_position_biotype+': '+str(round(value*100,2))+'% - ' 

				if type_gene_level_biotype == '':
					type_gene_level_biotype = 'missed peptide -'
				to_add = [type_peptide, peptide, type_gene_level_biotype[:-2], type_transcript_level_biotype[:-2], type_genomic_position_biotype[:-2], total_count_rna, total_count_ribo]
				data_final.append(to_add)


			# Peptide type
			count_gene_level_biotypes = Counter(gene_level_biotypes_type_peptide_rna)
			count_transcript_level_biotypes = Counter(transcript_level_biotypes_type_peptide_rna)
			count_genomic_position_biotypes = Counter(genomic_position_biotypes_type_peptide_rna)
			
			total_gene_level_biotype = len(gene_level_biotypes_type_peptide_rna)
			
			for gene_level_biotype, v in count_gene_level_biotypes.items():
				value =  v / (total_gene_level_biotype * 1.0)

				try:
					biotypes_peptides_rna_all_peptides['gene_level_biotype'][gene_level_biotype] += value
				except KeyError:
					biotypes_peptides_rna_all_peptides['gene_level_biotype'][gene_level_biotype] = value
				

			total_transcript_level_biotype = len(transcript_level_biotypes_type_peptide_rna)
			
			for transcript_level_biotype, v in count_transcript_level_biotypes.items():
				value =  v / (total_transcript_level_biotype*1.0)

				try:
					biotypes_peptides_rna_all_peptides['transcript_level_biotype'][transcript_level_biotype] += value
				except KeyError:
					biotypes_peptides_rna_all_peptides['transcript_level_biotype'][transcript_level_biotype] = value


			total_genomic_position_biotype = len(genomic_position_biotypes_type_peptide_rna)
			
			for genomic_position_biotype, v in count_genomic_position_biotypes.items():
				value =  v / (total_genomic_position_biotype*1.0)

				try:
					biotypes_peptides_rna_all_peptides['genomic_position_biotype'][genomic_position_biotype] += value
				except KeyError:
					biotypes_peptides_rna_all_peptides['genomic_position_biotype'][genomic_position_biotype] = value

			count_gene_level_biotypes = Counter(gene_level_biotypes_type_peptide_ribo)
			count_transcript_level_biotypes = Counter(transcript_level_biotypes_type_peptide_ribo)
			count_genomic_position_biotypes = Counter(genomic_position_biotypes_type_peptide_ribo)
			
			total_gene_level_biotype = len(gene_level_biotypes_type_peptide_ribo)
			
			for gene_level_biotype, v in count_gene_level_biotypes.items():
				value =  v / (total_gene_level_biotype * 1.0)

				try:
					biotypes_peptides_ribo_all_peptides['gene_level_biotype'][gene_level_biotype] += value
				except KeyError:
					biotypes_peptides_ribo_all_peptides['gene_level_biotype'][gene_level_biotype] = value
				

			total_transcript_level_biotype = len(transcript_level_biotypes_type_peptide_ribo)
			
			for transcript_level_biotype, v in count_transcript_level_biotypes.items():
				value =  v / (total_transcript_level_biotype*1.0)

				try:
					biotypes_peptides_ribo_all_peptides['transcript_level_biotype'][transcript_level_biotype] += value
				except KeyError:
					biotypes_peptides_ribo_all_peptides['transcript_level_biotype'][transcript_level_biotype] = value


			total_genomic_position_biotype = len(genomic_position_biotypes_type_peptide_ribo)
			
			for genomic_position_biotype, v in count_genomic_position_biotypes.items():
				value =  v / (total_genomic_position_biotype*1.0)

				try:
					biotypes_peptides_ribo_all_peptides['genomic_position_biotype'][genomic_position_biotype] += value
				except KeyError:
					biotypes_peptides_ribo_all_peptides['genomic_position_biotype'][genomic_position_biotype] = value
				


		# https://xlsxwriter.readthedocs.io/example_pandas_multiple.html
		self.df = pd.DataFrame(data, columns=['Peptide Type', 'Peptide','Alignment', 'Strand', 'Transcript', 'gene_level_biotype', 'transcript_level_biotype', 'genomic_position_biotype', 'reads count RNA', 'reads count Ribo'])
		self.df2 = pd.DataFrame(data_resume, columns=['Peptide Type', 'Peptide','Alignment', 'Strand', 'gene_level_biotype', 'transcript_level_biotype', 'genomic_position_biotype', 'reads count RNA', 'reads count Ribo'])
		self.df3 = pd.DataFrame(data_final, columns=['Peptide Type', 'Peptide', 'gene_level_biotype', 'transcript_level_biotype', 'genomic_position_biotype', 'Total reads count RNA', 'Total reads count Ribo'])

		self.write_xls_with_all_info_biotypes()

		self.draw_biotypes(biotypes_peptides_rna, self.path_to_output_folder+'plots/biotypes/genome_annotation/biotype_by_peptide_type_reads/by_level_biotype/', True)
		self.draw_biotypes(biotypes_peptides_ribo, self.path_to_output_folder+'plots/biotypes/genome_annotation/biotype_by_peptide_type_reads/by_level_biotype/', False)
		self.draw_biotypes(biotypes_peptides_global, self.path_to_output_folder+'plots/biotypes/genome_annotation/global_biotypes/by_peptide_type/', False, True)
		
		self.draw_biotypes_all_peptides(biotypes_peptides_all_peptides, self.path_to_output_folder+'plots/biotypes/genome_annotation/global_biotypes/all_peptides/', False, True)
		self.draw_biotypes_all_peptides(biotypes_peptides_rna_all_peptides,self.path_to_output_folder+'plots/biotypes/genome_annotation/biotype_by_peptide_type_reads/all_peptides/', True)
		self.draw_biotypes_all_peptides(biotypes_peptides_ribo_all_peptides, self.path_to_output_folder+'plots/biotypes/genome_annotation/biotype_by_peptide_type_reads/all_peptides/', False)
		
		if self.mode == 'translation':
			self.draw_correlation(data)
			

	def draw_correlation(self, data):

		#[type_peptide, peptide, alignment, strand, 'No Annotation', gene_level_biotype, transcript_level_biotype, genomic_position_biotype, count_rna, count_ribo]
		counts_ribo = []
		counts_rna = []
		peptide_type = []
		peptides = []

		for alignment in data:
			count_rna = math.log(alignment[-2]+1,10)
			count_ribo =  math.log(alignment[-1]+1,10)
			if count_ribo+count_rna  > 0:
				counts_ribo.append(count_ribo)
				counts_rna.append(count_rna)
				peptide_type.append(alignment[0])
				peptides.append(alignment[1])

		data = {'Type Peptide': peptide_type, 
				'Peptides': peptides,
        		'Read count Ribo': counts_ribo,
        		'Read count RNA': counts_rna} 
        		
		df_correlation = pd.DataFrame(data)
		plots.correlation(self.path_to_output_folder, self.name_exp, df_correlation)


	def draw_biotypes(self, biotypes_peptides, path_to_save, rna, global_ = False):

		others_non_coding = ['retained_intron', 'bidirectional_promoter_lncRNA', 'transcribed_unitary_pseudogene', 'transcribed_unprocessed_pseudogene', 'sense_overlapping','processed_pseudogene', 'unprocessed_pseudogene']
		others_protein_coding = ['IG_V_gene', 'TEC']

		organisation_labels = {'Protein-coding genes':['5UTR', '3UTR', 'In_frame', 'Frameshift', 'protein_coding', 'Junctions', 'Other coding regions'], 
							'Non-coding genes':['processed_transcript', 'nonsense_mediated_decay', 'antisense', 'Exons', 'lincRNA', 'Other non-coding regions'], 
							'Protein-coding transcripts':['5UTR', '3UTR', 'In_frame', 'Frameshift', 'protein_coding', 'CDS', 'Junctions', 'Other coding regions'], 
							'Non-coding transcripts':['processed_transcript', 'nonsense_mediated_decay', 'antisense', 'Exons', 'lincRNA', 'Other non-coding regions'], 
							'Intergenic Region':['Intergenic'], 
							'Intronic Region':['Introns']}

		for type_peptide, level_biotypes in biotypes_peptides.items():

			for level_biotype, info_level_biotype in level_biotypes.items():

				labels_in_type_peptide = {'Protein-coding genes':{}, 'Non-coding genes': {}, 'Protein-coding transcripts':{}, 'Non-coding transcripts': {}, 'Intergenic Region':{}, 'Intronic Region':{}}
				outer_labels = []
				outer_sizes = []
				intra_labels = []
				intra_sizes = []
				
				if len(info_level_biotype.values()) > 0:
					title = level_biotype+' '+type_peptide+' peptides_'

					for biotype, total_biotype in info_level_biotype.items():
						in_ = False
						for type_, types in organisation_labels.items():
							if biotype in types:
								if level_biotype == 'transcript_level_biotype' or level_biotype =='genomic_position_biotype':
									type_ = type_.replace('genes', 'transcripts')
								labels_in_type_peptide[type_][biotype] = total_biotype 
								in_ = True
								break
						if not in_:
							
							if biotype in others_non_coding:
								type_ = 'Non-coding genes'
								if level_biotype == 'transcript_level_biotype' or level_biotype == 'genomic_position_biotype':
									type_ = type_.replace('genes', 'transcripts')
								labels_in_type_peptide[type_]['Other non-coding regions'] = total_biotype 
							elif biotype in others_protein_coding:
								type_ = 'Protein-coding genes'
								if level_biotype == 'transcript_level_biotype' or level_biotype == 'genomic_position_biotype':
									type_ = type_.replace('genes', 'transcripts')
								labels_in_type_peptide[type_]['Other coding regions'] = total_biotype
							elif level_biotype =='genomic_position_biotype':
								labels_in_type_peptide['Protein-coding transcripts']['Junctions'] = total_biotype 
							else:
								print ('here ',biotype)

					for outer_label_type, intra_labels_type in organisation_labels.items():
						values_in = 0
						for intra_label in intra_labels_type:
							try:
								value = labels_in_type_peptide[outer_label_type][intra_label]
								intra_labels.append(intra_label)
								intra_sizes.append(value)
								values_in += value
							except:
								pass
						if values_in > 0:
							outer_labels.append(outer_label_type)
							outer_sizes.append(values_in)
					
					if not global_:
						if rna:
							plots.plot_pie(title+'rna', outer_labels, intra_labels, intra_sizes, outer_sizes,  path_to_save+level_biotype+'/', self.name_exp, type_peptide+'_rna_'+level_biotype)
						else:
							plots.plot_pie(title+'ribo', outer_labels, intra_labels, intra_sizes, outer_sizes,  path_to_save+level_biotype+'/', self.name_exp, type_peptide+'_ribo_'+level_biotype)
					else:
						plots.plot_pie(title+'rna', outer_labels, intra_labels, intra_sizes, outer_sizes,  path_to_save, self.name_exp, type_peptide+'_'+level_biotype)
						

	def draw_biotypes_all_peptides(self, biotypes_peptides, path_to_save, rna, global_=True):

		others_non_coding = ['retained_intron', 'bidirectional_promoter_lncRNA', 'transcribed_unitary_pseudogene', 'transcribed_unprocessed_pseudogene', 'sense_overlapping','processed_pseudogene', 'unprocessed_pseudogene']
		others_protein_coding = ['IG_V_gene', 'TEC']

		organisation_labels = {'Protein-coding genes':['5UTR', '3UTR', 'In_frame', 'Frameshift', 'protein_coding', 'Junctions', 'Other coding regions'], 
							'Non-coding genes':['processed_transcript', 'nonsense_mediated_decay', 'antisense', 'Exons', 'lincRNA', 'Other non-coding regions'], 
							'Protein-coding transcripts':['5UTR', '3UTR', 'In_frame', 'Frameshift', 'protein_coding', 'CDS', 'Junctions', 'Other coding regions'], 
							'Non-coding transcripts':['processed_transcript', 'nonsense_mediated_decay', 'antisense', 'Exons', 'lincRNA', 'Other non-coding regions'], 
							'Intergenic Region':['Intergenic'], 
							'Intronic Region':['Introns']}

		for level_biotype, info_level_biotype in biotypes_peptides.items():

			labels_in_type_peptide = {'Protein-coding genes':{}, 'Non-coding genes': {}, 'Protein-coding transcripts':{}, 'Non-coding transcripts': {}, 'Intergenic Region':{}, 'Intronic Region':{}}
			outer_labels = []
			outer_sizes = []
			intra_labels = []
			intra_sizes = []
			
			for biotype, total_biotype in info_level_biotype.items():
				title = level_biotype+' All peptides'
				in_ = False
				for type_, types in organisation_labels.items():
					if biotype in types:
						if level_biotype == 'transcript_level_biotype' or level_biotype =='genomic_position_biotype':
							type_ = type_.replace('genes', 'transcripts')
						labels_in_type_peptide[type_][biotype] = total_biotype 
						in_ = True
						break
				if not in_:
					
					if biotype in others_non_coding:
						type_ = 'Non-coding genes'
						if level_biotype == 'transcript_level_biotype' or level_biotype == 'genomic_position_biotype':
							type_ = type_.replace('genes', 'transcripts')
						labels_in_type_peptide[type_]['Other non-coding regions'] = total_biotype 
					elif biotype in others_protein_coding:
						type_ = 'Protein-coding genes'
						if level_biotype == 'transcript_level_biotype' or level_biotype == 'genomic_position_biotype':
							type_ = type_.replace('genes', 'transcripts')
						labels_in_type_peptide[type_]['Other coding regions'] = total_biotype
					elif level_biotype =='genomic_position_biotype':
						labels_in_type_peptide['Protein-coding transcripts']['Junctions'] = total_biotype 
					else:
						print ('here ',biotype)

			for outer_label_type, intra_labels_type in organisation_labels.items():
				values_in = 0
				for intra_label in intra_labels_type:
					try:
						value = labels_in_type_peptide[outer_label_type][intra_label]
						intra_labels.append(intra_label)
						intra_sizes.append(value)
						values_in += value
					except:
						pass
				if values_in > 0:
					outer_labels.append(outer_label_type)
					outer_sizes.append(values_in)
			
			if not global_:
				if rna:
					plots.plot_pie(title+'rna', outer_labels, intra_labels, intra_sizes, outer_sizes,  path_to_save, self.name_exp, 'All_peptides_rna_'+level_biotype)
				else:
					plots.plot_pie(title+'ribo', outer_labels, intra_labels, intra_sizes, outer_sizes,  path_to_save, self.name_exp, 'All_peptides_ribo_'+level_biotype)
			else:
				plots.plot_pie(title, outer_labels, intra_labels, intra_sizes, outer_sizes,  path_to_save, self.name_exp, 'All_peptides_'+level_biotype)



	def write_xls_with_all_info_biotypes(self):
		writer = pd.ExcelWriter(self.path_to_output_folder+'/res/Annotation_Biotypes.xlsx', engine='xlsxwriter')
		self.df.to_excel(writer, sheet_name='Transcript Annotation Biotypes')
		self.df2.to_excel(writer, sheet_name='Global Annotation Biotypes')
		self.df3.to_excel(writer, sheet_name='Peptide Final Annot. Biotypes')
		self.df4.to_excel(writer, sheet_name='ERE Annotation Biotypes')
		self.df5.to_excel(writer, sheet_name='Peptide final ERE Anno.Biotypes')
		writer.save()

		