import pickle
import pandas as pd
import utils.useful_functions as uf
import plotting.plots as plots
from collections import Counter
from genomics.biotype_genomic_information import BiotypeGenomicSearch


__author__ = "Maria Virginia Ruiz Cuevas"

class GetInformationBEDIntersection:

	def __init__(self, path_to_output_folder, mode, mouse):
		self.path_to_output_folder = path_to_output_folder

		if mode == 'translation':
			self.path_to_output_folder_bed_files = path_to_output_folder+'res_translation/BED_files/'
		else:
			self.path_to_output_folder_bed_files = self.path_to_output_folder+'res/BED_files/'

		path = self.path_to_output_folder +'/alignments/alignments_summary_information.pkl'
		self.alignments_summary_information = pd.read_pickle(path)
		self.mouse = mouse
		

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

					
		self.biotype_genomic_annotation_search = BiotypeGenomicSearch(self.peptides_intersected, genome_version, self.path_to_output_folder, self.mouse)
		self.information_final_biotypes_peptides = self.biotype_genomic_annotation_search.get_biotype_from_intersected_transcripts()
		
		with open(self.path_to_output_folder_bed_files+'/information_final_biotypes_peptides.dic', 'wb') as handle:
			pickle.dump(self.information_final_biotypes_peptides, handle, protocol=pickle.HIGHEST_PROTOCOL)
		

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


		with open(self.path_to_output_folder_bed_files+'/peptides_intersected_ere.dic', 'wb') as handle:
			pickle.dump(self.peptides_intersected_ere, handle, protocol=pickle.HIGHEST_PROTOCOL)


	def get_ribosome_profiling_transcripts_overlap(self, files_to_intersect):

		
		self.peptides_translated = {}
		self.transcripts_intersected = {}

		for sample, info_sample in files_to_intersect.items():
			bed_intersected = info_sample[2]
			
			dic = {}
			self.peptides_translated[sample] = dic
			# chr13	113208170	113208199	DEFGDSRRRW_chr13:113208170-113208199	-	chr13	StringTie	exon	113208024	113208650	1000	-	.	gene_id "STRG.85758"; transcript_id "STRG.85758.3"; exon_number "14"; reference_id "ENST00000375457.2"; ref_gene_id "ENSG00000126226.22"; ref_gene_name "PCID2"; cov "8.798036";	29
			# chr13	113208170	113208199	DEFGDSRRRW_chr13:113208170-113208199	-	chr13	StringTie	transcript	113177611	113208715	1000	-	.	gene_id "STRG.85758"; transcript_id "STRG.85758.1"; reference_id "ENST00000375479.6"; ref_gene_id "ENSG00000126226.22"; ref_gene_name "PCID2"; cov "43.302032"; FPKM "15.555397"; TPM "15.598426";	29
			with open(bed_intersected) as f:

				for index, line in enumerate(f):
					splitLine = line.strip().split('\t')
					chr = splitLine[0]
					start = splitLine[1]
					end = splitLine[2]
					peptide = splitLine[3].split('_')[0]
					position = splitLine[3].split('_')[1]
					strand_peptide = splitLine[4]
					strand_transcript = splitLine[11]
					

					key_peptide_position = chr+':'+start+'-'+end

					overlap = int(splitLine[-1])
					len_pos = int(end) - int(start)
					transcript = line.split(' transcript_id ')[1].split('\"')[1] 
					
					if 'transcript' in splitLine[7]:
						tpm = float(line.split(' TPM ')[1].split('\"')[1])
						self.transcripts_intersected[transcript+'_'+sample] = tpm

					key_peptide = splitLine[3]
					if overlap == len_pos and strand_peptide == strand_transcript and splitLine[7] == 'exon':
						split_key = key_peptide.split('_')[1].split('|')
						try:
							self.peptides_translated[sample][key_peptide][key_peptide_position].add(transcript)
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
							self.peptides_translated[sample][key_peptide] = dic
		
		with open(self.path_to_output_folder_bed_files+'/information_from_bed_files_intersected.dic', 'wb') as handle:
			pickle.dump(self.peptides_translated, handle, protocol=pickle.HIGHEST_PROTOCOL)

		with open(self.path_to_output_folder_bed_files+'/transcripts_intersected.dic', 'wb') as handle:
			pickle.dump(self.transcripts_intersected, handle, protocol=pickle.HIGHEST_PROTOCOL)
		
