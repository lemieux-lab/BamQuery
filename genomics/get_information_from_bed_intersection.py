import os, logging, time, pickle, multiprocessing, _thread, csv, math, copy, pysam
import pandas as pd
from pathos.multiprocessing import ProcessPool
import utils.useful_functions as uf
import plotting.plots as plots
from collections import Counter
from genomics.biotype_genomic_information import BiotypeGenomicSearch


path_to_lib = '/'.join(os.path.abspath(__file__).split('/')[:-3])+'/lib/'
__author__ = "Maria Virginia Ruiz Cuevas"

class GetInformationBEDIntersection:

	def __init__(self, path_to_output_folder):
		self.path_to_output_folder = path_to_output_folder
		self.path_to_output_folder_bed_files = self.path_to_output_folder+'res/BED_files/'
		path = self.path_to_output_folder +'/alignments/alignments_summary_information.pkl'
		self.alignments_summary_information = pd.read_pickle(path)
		

	def get_information_genomic_annotation(self, genome_version):

		file_to_read = self.path_to_output_folder_bed_files+'intersection_with_annotated_transcripts.bed'

		self.peptides_intersected = {}

		# chr13	99970439	99970465	AAAAPRPAL_chr13:99970439-99970465	+	chr13	HAVANA	transcript	99962964	99971909	.	+	.	gene_id "ENSG00000139800.8"; transcript_id "ENST00000267294.4"; gene_type "protein_coding"; gene_name "ZIC5"; transcript_type "protein_coding"; transcript_name "ZIC5-001"; level 1; protein_id "ENSP00000267294.3"; transcript_support_level "1"; tag "basic"; tag "appris_principal_1"; tag "exp_conf"; tag "CCDS"; ccdsid "CCDS9494.2"; havana_gene "OTTHUMG00000017280.3"; havana_transcript "OTTHUMT00000045623.3";	26
		with open(file_to_read) as f:

			for index, line in enumerate(f):
				splitLine = line.strip().split('\t')
				chr = splitLine[0]
				start = splitLine[1]
				end = splitLine[2]
				peptide = splitLine[3].split('_')[0]
				position = splitLine[3].split('_')[1]

				mcss = self.alignments_summary_information[(self.alignments_summary_information['Peptide'] == peptide) & (self.alignments_summary_information['Alignment'] == position)]['MCS'].values

				for mcs in mcss:
					key_peptide = splitLine[3]+'_'+mcs
					key_peptide_position = chr+':'+start+'-'+end

					strand_peptide = splitLine[4]
					strand_transcript = splitLine[11]
					
					if strand_peptide == strand_transcript :
						split_key = key_peptide.split('_')[1].split('|')
						transcript = ''
						transcript_support_level = ''
						
						intersection_number = int(splitLine[-1])

						try:
							transcript = line.split(' transcript_id ')[1].split('\"')[1]
						except IndexError:
							pass

						try:
							self.peptides_intersected[key_peptide][key_peptide_position].add(transcript)
						except KeyError:
							dic = {}
							transctipt_set = set()
							for key in split_key:
								if 'chr' not in key:
									key = chr+':'+key
								if key == key_peptide_position:
									transctipt_set.add(transcript)
									dic[key] = transctipt_set
								else:
									dic[key] = transctipt_set
							self.peptides_intersected[key_peptide] = dic

					
		self.biotype_genomic_annotation_search = BiotypeGenomicSearch(self.peptides_intersected, genome_version, self.path_to_output_folder)
		self.information_final_biotypes_peptides = self.biotype_genomic_annotation_search.get_biotype_from_intersected_transcripts()
		
		#with open(self.path_to_output_folder_bed_files+'/information_final_biotypes_peptides.dic', 'wb') as handle:
		#	pickle.dump(self.information_final_biotypes_peptides, handle, protocol=pickle.HIGHEST_PROTOCOL)
		

	def get_information_ERE_annotation(self):

		file_to_read = self.path_to_output_folder_bed_files+'intersection_with_annotated_EREs.bed'

		self.peptides_intersected_ere = {}
		# chr13	99970439	99970465	AAAAPRPAL_chr13:99970439-99970465	+	chr13	99970407	99970449	(GGC)n	41	+	10
		
		info_split_peptides = {}

		with open(file_to_read) as f:

			for index, line in enumerate(f):
				splitLine = line.strip().split('\t')
				chr = splitLine[0]
				start = splitLine[1]
				end = splitLine[2]
				peptide = splitLine[3].split('_')[0]
				number_overlap = int(splitLine[11])
				position = splitLine[3].split('_')[1]
				
				mcss = self.alignments_summary_information[(self.alignments_summary_information['Peptide'] == peptide) & (self.alignments_summary_information['Alignment'] == position)]['MCS'].values

				for mcs in mcss:
					key_peptide = splitLine[3]+'_'+mcs
					
					strand_peptide = splitLine[4]
					strand_transcript = splitLine[10]
					
					if number_overlap == ((len(peptide)*3) - 1):
						if strand_peptide == strand_transcript :
							key_position =  chr+':'+start+'-'+end+'_'+strand_peptide
							repName = splitLine[8]
								
							if '|' not in key_peptide: #((len(peptide)*3)/2.0):
								intersection_number = int(splitLine[11])
								
								try:
									peptide_info = self.peptides_intersected_ere[peptide]
									try:
										repName_set = peptide_info[key_peptide]
										repName_set.add(repName)
									except KeyError:
										repName_set = set()
										repName_set.add(repName)
										peptide_info[key_peptide] = repName_set
								except KeyError:
									repName_set = set()
									repName_set.add(repName)
									self.peptides_intersected_ere[peptide] = {key_peptide : repName_set}
							else:
								split_key = key_peptide.split('_')[1].split('|')

								for key in split_key:
									if 'chr' not in key:
										key = chr+':'+key+'_'+strand_peptide
									if key == key_position:
										break
								try:
									info_split_peptides[key_peptide][key][1] += number_overlap
								except KeyError:
									dic  = {}
									dic[key] = [repName, number_overlap]
									info_split_peptides[key_peptide] = dic
						else:
							key_position =  chr+':'+start+'-'+end+'_'+strand_peptide
							repName = splitLine[8]
								
							if '|' not in key_peptide: #((len(peptide)*3)/2.0):
								intersection_number = int(splitLine[11])
								
								try:
									peptide_info = self.peptides_intersected_ere[peptide]
									try:
										repName_set = peptide_info[key_peptide]
										repName_set.add('antisense_'+repName)
									except KeyError:
										repName_set = set()
										repName_set.add('antisense_'+repName)
										peptide_info[key_peptide] = repName_set
								except KeyError:
									repName_set = set()
									repName_set.add('antisense_'+repName)
									self.peptides_intersected_ere[peptide] = {key_peptide : repName_set}
							else:
								split_key = key_peptide.split('_')[1].split('|')

								for key in split_key:
									if 'chr' not in key:
										key = chr+':'+key+'_'+strand_peptide
									if key == key_position:
										break
								try:
									info_split_peptides[key_peptide][key][1] += number_overlap
								except KeyError:
									dic  = {}
									dic[key] = ['antisense_'+repName, number_overlap]
									info_split_peptides[key_peptide] = dic

					for key_peptide, info_split_peptide in info_split_peptides.items():
						
						peptide = key_peptide.split('_')[0]
						total_overlap = sum([ value[1] for key, value in info_split_peptide.items()])
						rep_names = set([ value[0] for key, value in info_split_peptide.items()])

						if total_overlap == ((len(peptide)*3) - 1):#total_overlap > ((len(peptide)*3)/2.0):
							try:
								peptide_info = self.peptides_intersected_ere[peptide]
								try:
									repName_set = peptide_info[key_peptide]
									repName_set = repName_set.union(rep_names)
								except KeyError:
									peptide_info[key_peptide] = rep_names
							except KeyError:
								self.peptides_intersected_ere[peptide] = {key_peptide : rep_names}


		#with open(self.path_to_output_folder_bed_files+'/peptides_intersected_ere.dic', 'wb') as handle:
		#	pickle.dump(self.peptides_intersected_ere, handle, protocol=pickle.HIGHEST_PROTOCOL)

		
