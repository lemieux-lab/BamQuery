import os, logging, time, pickle, multiprocessing, _thread, csv, math, copy, pysam
import pandas as pd
from pathos.multiprocessing import ProcessPool
import utils.useful_functions as uf
import plotting.plots as plots
from collections import Counter, OrderedDict
from genomics.get_information_from_bed_intersection import GetInformationBEDIntersection

path_to_lib = '/'.join(os.path.abspath(__file__).split('/')[:-3])+'/lib/'

__author__ = "Maria Virginia Ruiz Cuevas"

class BiotypeAssignation:

	def __init__(self, path_to_output_folder, name_exp, mode):
		self.path_to_output_folder = path_to_output_folder
		self.name_exp = name_exp
		self.mode = mode
		self.get_info_bed_files = GetInformationBEDIntersection(path_to_output_folder)
		self.get_info_bed_files.get_information_genomic_annotation()
		self.get_info_bed_files.get_information_ERE_annotation()
		with open(path_to_lib+'ERE_info.dic', 'rb') as handle:
			self.ere_info = pickle.load(handle)
		

	def get_biotypes(self, info_peptide_alignments, peptides_by_type, bam_files_list_rna, bam_files_list_ribo):
		
		data_ere = []
		data_gen = []
		data_gen_ere = []
		self.bam_files_list_rna = bam_files_list_rna
		self.bam_files_list_ribo = bam_files_list_ribo

		for type_peptide, peptides in peptides_by_type.items():

			for peptide in peptides:

				try:
					info_alignments_peptide = info_peptide_alignments[peptide]

					alignments = info_alignments_peptide[0]
					total_count_rna = info_alignments_peptide[1]
					total_count_ribo = info_alignments_peptide[2]
					
					# ('chr13:99970439-99970465', 'GCGGCGGCGGCACCCCGGCCAGCTCTT', '+', [0, 0], [])
					# information_final_biotypes_peptides[key_peptide] = {transcript: [gene_type, transcript_type, transcript_level_biotype]}
					for key_peptide_alignment, info_alignment in alignments.items():
						peptide = key_peptide_alignment.split('_')[0]
						position = key_peptide_alignment.split('_')[1]
						MCS = key_peptide_alignment.split('_')[2]
						strand = info_alignment[0]
						rna_bam_files = info_alignment[1]
						ribo_bam_files = info_alignment[2]
						
						to_add_ere = []
						to_add_gen = []
						to_add_gen_ere = []

						repNames = []
						repClass_s = []
						repFamilies = []

						key_aux = peptide+'_'+position

						try:
							rep_names = self.get_info_bed_files.peptides_intersected_ere[peptide][key_aux]
							
							for repName in rep_names:
								repClass = self.ere_info[repName][0]
								repFamily = self.ere_info[repName][1]
								repNames.append(repName)
								repClass_s.append(repClass)
								repFamilies.append(repFamily)

								to_add_ere = [type_peptide, peptide, position, MCS, strand, repName, repClass, repFamily]
								to_add_ere.extend(rna_bam_files)
								to_add_ere.extend(ribo_bam_files) 
								to_add_ere.extend([sum(rna_bam_files), sum(ribo_bam_files)])
								data_ere.append(to_add_ere)

						except KeyError:
							pass

						try:
							transcripts_biotypes = self.get_info_bed_files.information_final_biotypes_peptides[peptide][key_aux]

							for transcript, peptide_genomic_biotypes in transcripts_biotypes.items():
								gene_level_biotype = peptide_genomic_biotypes[0]
								transcript_level_biotype = peptide_genomic_biotypes[1]
								genomic_position_biotype = peptide_genomic_biotypes[2]
								if genomic_position_biotype == 'Exons':
									genomic_position_biotype = 'Non-coding Exons'

								to_add_gen = [type_peptide, peptide, position, MCS, strand, transcript, gene_level_biotype, transcript_level_biotype, genomic_position_biotype]
								to_add_gen_ere = [type_peptide, peptide, position, MCS, strand, transcript, gene_level_biotype, transcript_level_biotype, genomic_position_biotype, ",".join(repNames), ",".join(repClass_s), ",".join(repFamilies)]
								
								to_add_gen.extend(rna_bam_files)
								to_add_gen.extend(ribo_bam_files) 
								to_add_gen.extend([sum(rna_bam_files), sum(ribo_bam_files)])
								data_gen.append(to_add_gen)

								to_add_gen_ere.extend(rna_bam_files)
								to_add_gen_ere.extend(ribo_bam_files) 
								to_add_gen_ere.extend([sum(rna_bam_files), sum(ribo_bam_files)])
								data_gen_ere.append(to_add_gen_ere)

						except KeyError:
							transcript = 'No Annotation'
							gene_level_biotype = 'Intergenic'
							transcript_level_biotype = 'Intergenic'
							genomic_position_biotype = 'Intergenic'

							to_add_gen = [type_peptide, peptide, position, MCS, strand, transcript, gene_level_biotype, transcript_level_biotype, genomic_position_biotype]
							
							if len(repNames) > 0 :
								gene_level_biotype = ''
								transcript_level_biotype = ''
								genomic_position_biotype = ''
							to_add_gen_ere = [type_peptide, peptide, position, MCS, strand, transcript, gene_level_biotype, transcript_level_biotype, genomic_position_biotype, ",".join(repNames), ",".join(repClass_s), ",".join(repFamilies)]
							
							to_add_gen.extend(rna_bam_files)
							to_add_gen.extend(ribo_bam_files) 
							to_add_gen.extend([sum(rna_bam_files), sum(ribo_bam_files)])
							data_gen.append(to_add_gen)

							to_add_gen_ere.extend(rna_bam_files)
							to_add_gen_ere.extend(ribo_bam_files) 
							to_add_gen_ere.extend([sum(rna_bam_files), sum(ribo_bam_files)])
							data_gen_ere.append(to_add_gen_ere)

				except KeyError:
					info_alignments_peptide = []

		columns_ere = ['Peptide Type', 'Peptide','Alignment', 'MCS', 'Strand', 'ERE name', 'ERE class', 'ERE family']
		columns_ere.extend(list(self.bam_files_list_rna.keys()))
		columns_ere.extend(list(self.bam_files_list_ribo.keys()))
		columns_ere.extend(['Total reads count RNA', 'Total reads count Ribo'])

		columns_gen = ['Peptide Type', 'Peptide','Alignment', 'MCS', 'Strand', 'Transcript', 'gene_level_biotype', 'transcript_level_biotype', 'genomic_position_biotype']
		columns_gen.extend(list(self.bam_files_list_rna.keys()))
		columns_gen.extend(list(self.bam_files_list_ribo.keys()))
		columns_gen.extend(['Total reads count RNA', 'Total reads count Ribo'])

		columns_gen_ere = ['Peptide Type', 'Peptide','Alignment', 'MCS', 'Strand', 'Transcript', 'gene_level_biotype', 'transcript_level_biotype', 'genomic_position_biotype', 'ERE name', 'ERE class', 'ERE family']
		columns_gen_ere.extend(list(self.bam_files_list_rna.keys()))
		columns_gen_ere.extend(list(self.bam_files_list_ribo.keys()))
		columns_gen_ere.extend(['Total reads count RNA', 'Total reads count Ribo'])

		self.df1 = pd.DataFrame(data_ere, columns = columns_ere)
		self.df2 = pd.DataFrame(data_gen, columns = columns_gen)
		self.df3 = pd.DataFrame(data_gen_ere, columns = columns_gen_ere)
		self.write_xls_with_all_info_biotypes()
		

	def get_global_annotation(self):

		# ERE Full
		groupby_columns = ['Peptide Type', 'Peptide','Alignment', 'Strand', 'ERE name', 'ERE class', 'ERE family']
		bam_files_columns = list(self.bam_files_list_rna.keys())
		bam_files_columns.extend(list(self.bam_files_list_ribo.keys()))
		bam_files_columns.extend(['Total reads count RNA', 'Total reads count Ribo'])
		
		self.df_total_by_position_ere = self.df1.groupby(groupby_columns)[bam_files_columns].sum().reset_index()
		self.df_total_by_position_ere = pd.DataFrame(self.df_total_by_position_ere, columns = groupby_columns+bam_files_columns)
		
		# Genomic Annotation Full
		groupby_columns = ['Peptide Type', 'Peptide','Alignment', 'Strand', 'Transcript', 'gene_level_biotype', 'transcript_level_biotype', 'genomic_position_biotype']
		self.df_total_by_position_gen = self.df2.groupby(groupby_columns)[bam_files_columns].sum().reset_index()
		self.df_total_by_position_gen = pd.DataFrame(self.df_total_by_position_gen, columns = groupby_columns+bam_files_columns)
		
		groupby_columns = ['Peptide Type', 'Peptide','Alignment', 'Strand', 'Transcript', 'gene_level_biotype', 'transcript_level_biotype', 'genomic_position_biotype', 'ERE name', 'ERE class', 'ERE family']
		self.df_total_by_position_gen_ere = self.df3.groupby(groupby_columns)[bam_files_columns].sum().reset_index()
		self.df_total_by_position_gen_ere = pd.DataFrame(self.df_total_by_position_gen_ere, columns = groupby_columns+bam_files_columns)
		
		self.get_global_annotation_gen(bam_files_columns)
		self.get_global_annotation_ere(bam_files_columns)
		self.get_final_annotation_peptide_gen(bam_files_columns)
		self.get_final_annotation_peptide_ere(bam_files_columns)

		self.get_genomic_and_ere_annotation(bam_files_columns)
		self.write_xls_with_consensus_biotypes()

	def get_global_annotation_gen(self, bam_files_columns):

		groupby_columns = ['Peptide Type', 'Peptide', 'Alignment', 'Strand', 'gene_level_biotype']
		df_global_gen_gene = self.df_total_by_position_gen.groupby(groupby_columns)['gene_level_biotype'].size()

		groupby_columns = ['Peptide Type', 'Peptide', 'Alignment', 'Strand', 'transcript_level_biotype']
		df_global_gen_transcript = self.df_total_by_position_gen.groupby(groupby_columns)['transcript_level_biotype'].size()

		groupby_columns = ['Peptide Type', 'Peptide', 'Alignment', 'Strand', 'genomic_position_biotype']
		df_global_gen_pos = self.df_total_by_position_gen.groupby(groupby_columns)['genomic_position_biotype'].size()
		
		groupby_columns = ['Peptide Type', 'Peptide','Alignment', 'Strand']
		df_global_gen_construction = self.df_total_by_position_gen.groupby(groupby_columns+bam_files_columns)
		
		group = [df_global_gen_gene, df_global_gen_transcript, df_global_gen_pos]
		data_global_annotation_gen = self.get_consensus(df_global_gen_construction, group, len(groupby_columns))

		groupby_columns = ['Peptide Type', 'Peptide','Alignment', 'Strand', 'gene_level_biotype', 'transcript_level_biotype', 'genomic_position_biotype']
		self.df_global_gen = pd.DataFrame(data_global_annotation_gen, columns = groupby_columns+bam_files_columns)

	def get_global_annotation_ere(self, bam_files_columns):

		groupby_columns = ['Peptide Type', 'Peptide','Alignment', 'Strand', 'ERE class']
		df_global_ere_class = self.df_total_by_position_ere.groupby(groupby_columns)['ERE class'].size()

		groupby_columns = ['Peptide Type', 'Peptide', 'Alignment', 'Strand', 'ERE family']
		df_global_ere_family = self.df_total_by_position_ere.groupby(groupby_columns)['ERE family'].size()

		groupby_columns = ['Peptide Type', 'Peptide','Alignment', 'Strand']
		df_global_ere_construction = self.df_total_by_position_ere.groupby(groupby_columns+bam_files_columns)
		
		group = [df_global_ere_class, df_global_ere_family]
		data_global_annotation_ere = self.get_consensus(df_global_ere_construction, group, len(groupby_columns))

		groupby_columns = ['Peptide Type', 'Peptide','Alignment', 'Strand', 'ERE class', 'ERE family']
		self.df_global_ere = pd.DataFrame(data_global_annotation_ere, columns = groupby_columns+bam_files_columns)

	def get_final_annotation_peptide_gen(self, bam_files_columns):

		groupby_columns = ['Peptide Type', 'Peptide', 'gene_level_biotype']
		df_global_gen_gene = self.df_total_by_position_gen.groupby(groupby_columns)['gene_level_biotype'].size()

		groupby_columns = ['Peptide Type', 'Peptide', 'transcript_level_biotype']
		df_global_gen_transcript = self.df_total_by_position_gen.groupby(groupby_columns)['transcript_level_biotype'].size()

		groupby_columns = ['Peptide Type', 'Peptide', 'genomic_position_biotype']
		df_global_gen_pos = self.df_total_by_position_gen.groupby(groupby_columns)['genomic_position_biotype'].size()
		
		groupby_columns = ['Peptide Type', 'Peptide']
		df_final_gen_sum = self.df_global_gen.groupby(groupby_columns)[bam_files_columns].sum()
		df_final_gen_sum = df_final_gen_sum.groupby(groupby_columns)

		df_final_gen_construction = self.df_total_by_position_gen.groupby(groupby_columns)
		
		group = [df_global_gen_gene, df_global_gen_transcript, df_global_gen_pos]
		data_global_annotation_gen = self.get_final_consensus(df_final_gen_construction, df_final_gen_sum, group, len(groupby_columns))

		groupby_columns = ['Peptide Type', 'Peptide', 'gene_level_biotype', 'transcript_level_biotype', 'genomic_position_biotype']
		self.df_final_peptide_gen = pd.DataFrame(data_global_annotation_gen, columns = groupby_columns+bam_files_columns)

	def get_final_annotation_peptide_ere(self, bam_files_columns):

		groupby_columns = ['Peptide Type', 'Peptide', 'ERE class']
		df_global_ere_class = self.df_total_by_position_ere.groupby(groupby_columns)['ERE class'].size()

		groupby_columns = ['Peptide Type', 'Peptide', 'ERE family']
		df_global_ere_family = self.df_total_by_position_ere.groupby(groupby_columns)['ERE family'].size()

		groupby_columns = ['Peptide Type', 'Peptide']
		df_final_ere_sum = self.df_global_ere.groupby(groupby_columns)[bam_files_columns].sum()
		df_final_ere_sum = df_final_ere_sum.groupby(groupby_columns)
		
		df_final_ere_construction = self.df_total_by_position_ere.groupby(groupby_columns)
		
		group = [df_global_ere_class, df_global_ere_family]
		data_global_annotation_ere = self.get_final_consensus(df_final_ere_construction, df_final_ere_sum, group, len(groupby_columns))

		groupby_columns = ['Peptide Type', 'Peptide','ERE class', 'ERE family']
		self.df_final_peptide_ere = pd.DataFrame(data_global_annotation_ere, columns = groupby_columns+bam_files_columns)

	def get_consensus(self, df_global, group, key_length):

		df_global_annotation = []

		for key, items in df_global:

			to_add = []
			total_items = len(items)
			new_key = key[:key_length]
			to_add.extend(list(new_key))
			
			count_bamfiles = key[key_length:]
			consensus = ['']*len(group)

			for index, df in enumerate(group):

				res = df[new_key].sort_values(ascending=False)
				count_string = []

				for type_, count in res.items():

					percentage = str(round((count/(total_items*1.0)*100),2))
					count_string.append(type_+': '+percentage+'%')
				consensus[index] = ' - '.join(count_string)

			to_add.extend(consensus)
			to_add.extend(list(count_bamfiles))
			df_global_annotation.append(to_add)

		return df_global_annotation

	def get_final_consensus(self, df_global, df_final_sum, group, key_length):

		df_global_annotation = []

		for key, items in df_global:

			to_add = []
			total_items = len(items)
			new_key = key[:key_length]
			to_add.extend(list(new_key))
			
			count_bamfiles = df_final_sum.get_group(key).values.tolist()[0]
			consensus = ['']*len(group)
			
			for index, df in enumerate(group):

				res = df[new_key].sort_values(ascending=False)
				count_string = []

				for type_, count in res.items():

					percentage = str(round((count/(total_items*1.0)*100),2))
					count_string.append(type_+': '+percentage+'%')
				consensus[index] = ' - '.join(count_string)

			to_add.extend(consensus)
			to_add.extend(list(count_bamfiles))
			df_global_annotation.append(to_add)

		return df_global_annotation

	def get_genomic_and_ere_annotation(self, bam_files_columns):

		df_consensus_annotation = []
		groupby_columns = ['Peptide Type', 'Peptide', 'Alignment', 'Strand']
		group_by_gen_ere = self.df_total_by_position_gen_ere.groupby(groupby_columns)
		
		groupby_columns = ['Peptide Type', 'Peptide']
		df_final_gen_sum = self.df_global_gen.groupby(groupby_columns)[bam_files_columns].sum()
		df_final_gen_sum = df_final_gen_sum.groupby(groupby_columns)

		peptides_visited = {}

		biotypes_by_peptide_type = {}
		biotypes_all_peptides = {}
		counts_reads_rna_ribo = {}

		def insert_dic(biotype, sum_percentage):
			try:
				biotypes = peptides_visited[(peptide_type, peptide)][0]
				try:
					biotypes[biotype] += sum_percentage
				except KeyError:
					biotypes[biotype] = sum_percentage
			except KeyError:
				dic = {}
				dic[biotype] = sum_percentage
				peptides_visited[(peptide_type, peptide)] = [dic]

		for key, items in group_by_gen_ere:
			peptide_type = key[0]
			peptide = key[1]
			alignment = key[2]
			count_bamfiles = df_final_gen_sum.get_group((peptide_type, peptide)).values.tolist()[0]
			total_items_pep = 0
			
			type_gen = list(filter(lambda a:a != '', items['genomic_position_biotype'].values.tolist()))
			type_ere = list(filter(lambda a:a != '', items['ERE class'].values.tolist()))

			total_items = len(type_gen) + len(type_ere)
			sum_percentage = 1 / (total_items*1.0)
			
			for genomic_position_biotype in type_gen:
				if genomic_position_biotype != '' :
					insert_dic(genomic_position_biotype, sum_percentage)

			for ere_class in type_ere:
				if ere_class != '' :
					insert_dic(ere_class, sum_percentage)
			if len(peptides_visited[(peptide_type, peptide)]) == 1:
				peptides_visited[(peptide_type, peptide)].append(count_bamfiles)

		for key, information_peptide in peptides_visited.items():
			
			to_add = []
			peptide_type = key[0]
			peptide = key[1]
			total_alignments_pep = sum(peptides_visited[(peptide_type, peptide)][0].values())
			dic = peptides_visited[(peptide_type, peptide)][0]
			count_bamfiles = peptides_visited[(peptide_type, peptide)][1]
			
			count_rna = count_bamfiles[-2]
			count_ribo = count_bamfiles[-1]

			count_string = []

			for type_ in sorted(dic, key=dic.get, reverse=True):
				count = dic[type_]
				count = count/(total_alignments_pep*1.0)
				percentage = str(round(count*100,2))
				count_string.append(type_+': '+percentage+'%')
				
				try:
					types = biotypes_by_peptide_type[peptide_type]
					try:
						types[type_] += count
					except KeyError:
						types[type_] = count
				except KeyError:
					biotypes_by_peptide_type[peptide_type] = {type_ : count}

				try:
					biotypes_all_peptides[type_] += count
				except KeyError:
					biotypes_all_peptides[type_] = count

			try:
				counts_reads_rna_ribo[peptide_type][peptide] = [count_rna, count_ribo]
			except KeyError:
				dic = {}
				dic[peptide] = [count_rna, count_ribo]
				counts_reads_rna_ribo[peptide_type] = dic
			
			consensus = ' - '.join(count_string)
			to_add = [peptide_type, peptide, consensus]
			to_add.extend(count_bamfiles)
			df_consensus_annotation.append(to_add)

		groupby_columns = ['Peptide Type', 'Peptide', 'Consenssus']
		self.df_consensus_annotation = pd.DataFrame(df_consensus_annotation, columns = groupby_columns+bam_files_columns)

		self.draw_biotypes(biotypes_by_peptide_type, self.path_to_output_folder+'plots/biotypes/genome_and_ERE_annotation/by_peptide_type/', False)
		self.draw_biotypes(biotypes_all_peptides, self.path_to_output_folder+'plots/biotypes/genome_and_ERE_annotation/all_peptides/', True)
		
		if self.mode == 'translation':
			self.draw_correlation(counts_reads_rna_ribo)
			

	def draw_biotypes(self, biotypes_peptides, path_to_save, global_):

		others_non_coding = ['retained_intron', 'bidirectional_promoter_lncRNA', 'transcribed_unitary_pseudogene', 'transcribed_unprocessed_pseudogene', 'sense_overlapping','processed_pseudogene', 'unprocessed_pseudogene']
		others_protein_coding = ['IG_V_gene', 'TEC']
		others_ere = ['DNA','RC', 'RNA','Satellite','Simple_repeat','Unknown', 'Low_complexity', 'rRNA','scRNA','snRNA','srpRNA','tRNA']
		

		organisation_labels = {'Protein-coding Regions':['5UTR', '3UTR', 'In_frame', 'Frameshift', 'protein_coding', 'CDS', 'Junctions', 'Other coding regions'], 
							'Non-coding Regions':['processed_transcript', 'nonsense_mediated_decay', 'antisense', 'Non-coding Exons', 'lincRNA', 'Other non-coding regions'], 
							'Intergenic Regions':['Intergenic'], 
							'Intronic Regions':['Introns'],
							'EREs':['LINE', 'LTR','Retroposon','SINE', 'Other EREs']}

		def plot_biotype(biotypes, name):
			title = name+' peptides'
			labels_in_type_peptide = {'Protein-coding genes':{}, 'Non-coding genes': {}, 'Protein-coding Regions':{}, 'Non-coding Regions': {}, 'Protein-coding transcripts':{}, 'Non-coding transcripts': {}, 'Intergenic Regions':{}, 'Intronic Regions':{}, 'EREs':{}}

			
			for biotype, total_biotype in biotypes.items():

				outer_labels = []
				outer_sizes = []
				intra_labels = []
				intra_sizes = []
				
				in_ = False
				for type_, types in organisation_labels.items():
					if biotype in types:
						labels_in_type_peptide[type_][biotype] = total_biotype 
						in_ = True
						break

				if not in_:
						
					if biotype in others_non_coding:
						type_ = 'Non-coding Regions'
						labels_in_type_peptide[type_]['Other non-coding regions'] = total_biotype 

					elif biotype in others_protein_coding:
						type_ = 'Protein-coding Regions'
						labels_in_type_peptide[type_]['Other coding regions'] = total_biotype

					elif biotype in others_ere:
						type_ = 'EREs'
						labels_in_type_peptide[type_]['Other EREs'] = total_biotype

					else:
						labels_in_type_peptide['Protein-coding Regions']['Junctions'] = total_biotype 
						

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
					
				plots.plot_pie(title, outer_labels, intra_labels, intra_sizes, outer_sizes, path_to_save, self.name_exp, name)
		
		if not global_:
			
			for type_peptide, biotypes in biotypes_peptides.items():
				plot_biotype(biotypes, type_peptide)
		else:
			plot_biotype(biotypes_peptides, 'All_peptides')

			
	def prepare_info_to_draw_biotypes(self, info_peptide_alignments, peptides_by_type):

		biotypes_peptides_rna = {}
		biotypes_peptides_ribo = {}
		biotypes_peptides_rna_all_peptides = {'gene_level_biotype':{}, 'transcript_level_biotype':{}, 'genomic_position_biotype':{}}
		biotypes_peptides_ribo_all_peptides = {'gene_level_biotype':{}, 'transcript_level_biotype':{}, 'genomic_position_biotype':{}}
		biotypes_peptides_all_peptides = {'gene_level_biotype':{}, 'transcript_level_biotype':{}, 'genomic_position_biotype':{}}

		biotypes_unified_all_peptides = {'Genomic_and_ERE_Annotations_biotype':{}}
		biotypes_peptides_global = {}
		
		biotypes_peptides_global_unified = {}

		data = []
		data_resume = []
		data_final = []

		data_unified = []
		data_unified_final = []
		
		for type_peptide, peptides in peptides_by_type.items():
			biotypes_peptides_rna[type_peptide] = {'gene_level_biotype':{}, 'transcript_level_biotype':{}, 'genomic_position_biotype':{}}
			biotypes_peptides_ribo[type_peptide] = {'gene_level_biotype':{}, 'transcript_level_biotype':{}, 'genomic_position_biotype':{}}

			biotypes_peptides_global[type_peptide] = {'gene_level_biotype':{}, 'transcript_level_biotype':{}, 'genomic_position_biotype':{}}
			biotypes_peptides_global_unified[type_peptide] = {'Genomic_and_ERE_Annotations_biotype':{}}

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
				ERE_class_genomic_positions_biotypes_peptide = []

				gene_level_biotypes = []
				transcript_level_biotypes = []
				genomic_position_biotypes = []
				ERE_class_genomic_positions_biotypes = []

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

						key_peptide_position_ere = alignment+'_'+strand
						ERE_class_genomic_positions_biotypes = []

						total_count_rna += count_rna
						total_count_ribo += count_ribo
						
						try:
							transcripts_intersected = self.information_final_biotypes_peptides[key]
							try:
								repName = list(self.get_info_bed_files.peptides_intersected_ere[peptide][key_peptide_position_ere])
							except KeyError:
								repName = []

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
								ERE_class_genomic_positions_biotypes.append(genomic_position_biotype)

								to_add = [type_peptide, peptide, alignment, strand, position, gene_level_biotype, transcript_level_biotype, genomic_position_biotype, count_rna, count_ribo]
								data.append(to_add)

								if repName:
									for name in repName:
										repClass = self.ere_info[name][0]
										repFamily = self.ere_info[name][1]
										to_add_ere = [type_peptide, peptide, alignment, strand, position, gene_level_biotype, transcript_level_biotype, genomic_position_biotype, name, repFamily, repClass, count_rna, count_ribo]
										ERE_class_genomic_positions_biotypes.append(repClass)
										data_unified.append(to_add_ere)
								else:
									repClass = ''
									repFamily = ''
									name = ''
									to_add_ere = [type_peptide, peptide, alignment, strand, position, gene_level_biotype, transcript_level_biotype, genomic_position_biotype, name, repFamily, repClass, count_rna, count_ribo]
									data_unified.append(to_add_ere)
								
							count_gene_level_biotypes = Counter(gene_level_biotypes).most_common()
							count_transcript_level_biotypes = Counter(transcript_level_biotypes).most_common()
							count_genomic_position_biotypes = Counter(genomic_position_biotypes).most_common()
							count_ERE_class_genomic_postions_biotypes = Counter(ERE_class_genomic_positions_biotypes).most_common()

							total_gene_level_biotype = len(gene_level_biotypes)
							type_gene_level_biotype = ''
							
							for gene_level_biotype in count_gene_level_biotypes:
								
								v = gene_level_biotype[1]
								gene_level_biotype = gene_level_biotype[0]

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
							
							for transcript_level_biotype in count_transcript_level_biotypes:

								v = transcript_level_biotype[1]
								transcript_level_biotype = transcript_level_biotype[0]

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
							
							for genomic_position_biotype in count_genomic_position_biotypes:
								v = genomic_position_biotype[1]
								genomic_position_biotype = genomic_position_biotype[0]

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
							
							
							for ERE_class_genomic_position_biotype in count_ERE_class_genomic_postions_biotypes:
								v = ERE_class_genomic_position_biotype[1]
								ERE_class_genomic_position_biotype = ERE_class_genomic_position_biotype[0]

								ERE_class_genomic_positions_biotypes_peptide.append(ERE_class_genomic_position_biotype)

								try:
									biotypes_peptides_global_unified[type_peptide]['Genomic_and_ERE_Annotations_biotype'][ERE_class_genomic_position_biotype] += value
								except KeyError:
									biotypes_peptides_global_unified[type_peptide]['Genomic_and_ERE_Annotations_biotype'][ERE_class_genomic_position_biotype] = value


							to_add = [type_peptide, peptide, alignment, strand, type_gene_level_biotype[:-2], type_transcript_level_biotype[:-2], type_genomic_position_biotype[:-2], count_rna, count_ribo]
							
							data_resume.append(to_add)
	
						except KeyError:

							try:
								repName = list(self.peptides_intersected_ere[peptide][key_peptide_position_ere])
							except KeyError:
								repName = []

							biotype_type = 'Intergenic'
							gene_level_biotype = 'Intergenic'
							transcript_level_biotype = 'Intergenic'
							genomic_position_biotype = biotype_type

							gene_level_biotypes = [gene_level_biotype]
							transcript_level_biotypes = [transcript_level_biotype]
							genomic_position_biotypes = [genomic_position_biotype]

							gene_level_biotypes_peptide.append(gene_level_biotype)
							transcript_level_biotypes_peptide.append(transcript_level_biotype)
							genomic_position_biotypes_peptide.append(genomic_position_biotype)

							to_add = [type_peptide, peptide, alignment, strand, 'No Annotation', gene_level_biotype, transcript_level_biotype, genomic_position_biotype, count_rna, count_ribo]
							data.append(to_add)
							
							to_add = [type_peptide, peptide, alignment, strand, gene_level_biotype, transcript_level_biotype, genomic_position_biotype, count_rna, count_ribo]
							data_resume.append(to_add)

							if repName:
								gene_level_biotype_aux = ''
								transcript_level_biotype_aux = ''
								genomic_position_biotype_aux = ''
								position_aux = 'No Annotation'

								for name in repName:
									repClass = self.ere_info[name][0]
									repFamily = self.ere_info[name][1]
									to_add_ere = [type_peptide, peptide, alignment, strand, position_aux, gene_level_biotype_aux, transcript_level_biotype_aux, genomic_position_biotype_aux, name, repFamily, repClass, count_rna, count_ribo]
									ERE_class_genomic_positions_biotypes_peptide.append(repClass)
									ERE_class_genomic_positions_biotypes.append(repClass)
									data_unified.append(to_add_ere)

									try:
										biotypes_peptides_global_unified[type_peptide]['Genomic_and_ERE_Annotations_biotype'][repClass] += 1
									except KeyError:
										biotypes_peptides_global_unified[type_peptide]['Genomic_and_ERE_Annotations_biotype'][repClass] = 1
							else:
								gene_level_biotype_aux = biotype_type
								transcript_level_biotype_aux = biotype_type
								genomic_position_biotype_aux = biotype_type
								position_aux = 'No Annotation'

								try:
									biotypes_peptides_global_unified[type_peptide]['Genomic_and_ERE_Annotations_biotype'][biotype_type] += 1
								except KeyError:
									biotypes_peptides_global_unified[type_peptide]['Genomic_and_ERE_Annotations_biotype'][biotype_type] = 1

								repClass = ''
								repFamily = ''
								name = ''

								ERE_class_genomic_positions_biotypes_peptide.append(biotype_type)
								ERE_class_genomic_positions_biotypes.append(biotype_type)

								to_add_ere = [type_peptide, peptide, alignment, strand, position_aux, gene_level_biotype_aux, transcript_level_biotype_aux, genomic_position_biotype_aux, name, repFamily, repClass, count_rna, count_ribo]
								data_unified.append(to_add_ere)


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
				if len(alignments) == 1 :
					count_gene_level_biotypes = Counter(gene_level_biotypes).most_common()
					total_gene_level_biotype = len(gene_level_biotypes)

					count_transcript_level_biotypes = Counter(transcript_level_biotypes).most_common()
					total_transcript_level_biotype = len(transcript_level_biotypes)

					count_genomic_position_biotypes = Counter(genomic_position_biotypes).most_common()
					total_genomic_position_biotype = len(genomic_position_biotypes)

					count_ERE_class_genomic_postions_biotypes = Counter(ERE_class_genomic_positions_biotypes).most_common()
					total_ERE_class_genomic_positions_biotype = len(ERE_class_genomic_positions_biotypes)

				else:
					count_gene_level_biotypes = Counter(gene_level_biotypes_peptide).most_common()
					total_gene_level_biotype = len(gene_level_biotypes_peptide)

					count_transcript_level_biotypes = Counter(transcript_level_biotypes_peptide).most_common()
					total_transcript_level_biotype = len(transcript_level_biotypes_peptide)

					count_genomic_position_biotypes = Counter(genomic_position_biotypes_peptide).most_common()
					total_genomic_position_biotype = len(genomic_position_biotypes_peptide)

					count_ERE_class_genomic_postions_biotypes = Counter(ERE_class_genomic_positions_biotypes_peptide).most_common()
					total_ERE_class_genomic_positions_biotype = len(ERE_class_genomic_positions_biotypes_peptide)

				
				type_gene_level_biotype = ''
				
				for gene_level_biotype in count_gene_level_biotypes:
					v = gene_level_biotype[1]
					gene_level_biotype = gene_level_biotype[0]

					value =  v / (total_gene_level_biotype*1.0)
					try:
						biotypes_peptides_all_peptides['gene_level_biotype'][gene_level_biotype] += value
					except KeyError:
						biotypes_peptides_all_peptides['gene_level_biotype'][gene_level_biotype] = value

					type_gene_level_biotype += gene_level_biotype+': '+str(round(value*100,2))+'% - ' 
					
				type_transcript_level_biotype = ''
				
				for transcript_level_biotype in count_transcript_level_biotypes :
					v = transcript_level_biotype[1]
					transcript_level_biotype = transcript_level_biotype[0]

					value =  v / (total_transcript_level_biotype*1.0)
					try:
						biotypes_peptides_all_peptides['transcript_level_biotype'][transcript_level_biotype] += value
					except KeyError:
						biotypes_peptides_all_peptides['transcript_level_biotype'][transcript_level_biotype] = value

					type_transcript_level_biotype += transcript_level_biotype+': '+str(round(value*100,2))+'% - ' 
					
				
				type_genomic_position_biotype = ''
				
				for genomic_position_biotype in count_genomic_position_biotypes:
					v = genomic_position_biotype[1]
					genomic_position_biotype = genomic_position_biotype[0]

					value =  v / (total_genomic_position_biotype*1.0)
					try:
						biotypes_peptides_all_peptides['genomic_position_biotype'][genomic_position_biotype] += value
					except KeyError:
						biotypes_peptides_all_peptides['genomic_position_biotype'][genomic_position_biotype] = value

					type_genomic_position_biotype += genomic_position_biotype+': '+str(round(value*100,2))+'% - ' 

				
				type_ERE_class_genomic_positions_biotype = ''
				
				for ERE_class_genomic_positions_biotype in count_ERE_class_genomic_postions_biotypes:
					v = ERE_class_genomic_positions_biotype[1]
					ERE_class_genomic_positions_biotype = ERE_class_genomic_positions_biotype[0]

					value =  v / (total_ERE_class_genomic_positions_biotype*1.0)
					try:
						biotypes_unified_all_peptides['Genomic_and_ERE_Annotations_biotype'][ERE_class_genomic_positions_biotype] += value
					except KeyError:
						biotypes_unified_all_peptides['Genomic_and_ERE_Annotations_biotype'][ERE_class_genomic_positions_biotype] = value

					type_ERE_class_genomic_positions_biotype += ERE_class_genomic_positions_biotype+': '+str(round(value*100,2))+'% - ' 

				if type_gene_level_biotype == '':
					type_gene_level_biotype = 'missed peptide -'
				
				to_add = [type_peptide, peptide, type_gene_level_biotype[:-2], type_transcript_level_biotype[:-2], type_genomic_position_biotype[:-2], total_count_rna, total_count_ribo]
				data_final.append(to_add)

				to_add = [type_peptide, peptide, type_ERE_class_genomic_positions_biotype[:-2], total_count_rna, total_count_ribo]
				data_unified_final.append(to_add)


			# Peptide type
			count_gene_level_biotypes = Counter(gene_level_biotypes_type_peptide_rna).most_common()
			count_transcript_level_biotypes = Counter(transcript_level_biotypes_type_peptide_rna).most_common()
			count_genomic_position_biotypes = Counter(genomic_position_biotypes_type_peptide_rna).most_common()
			
			total_gene_level_biotype = len(gene_level_biotypes_type_peptide_rna)
			
			for gene_level_biotype in count_gene_level_biotypes:

				v = gene_level_biotype[1]
				gene_level_biotype = gene_level_biotype[0]

				value =  v / (total_gene_level_biotype * 1.0)

				try:
					biotypes_peptides_rna_all_peptides['gene_level_biotype'][gene_level_biotype] += value
				except KeyError:
					biotypes_peptides_rna_all_peptides['gene_level_biotype'][gene_level_biotype] = value
				

			total_transcript_level_biotype = len(transcript_level_biotypes_type_peptide_rna)
			
			for transcript_level_biotype in count_transcript_level_biotypes :

				v = transcript_level_biotype[1]
				transcript_level_biotype = transcript_level_biotype[0]
				value =  v / (total_transcript_level_biotype*1.0)

				try:
					biotypes_peptides_rna_all_peptides['transcript_level_biotype'][transcript_level_biotype] += value
				except KeyError:
					biotypes_peptides_rna_all_peptides['transcript_level_biotype'][transcript_level_biotype] = value


			total_genomic_position_biotype = len(genomic_position_biotypes_type_peptide_rna)
			
			for genomic_position_biotype  in count_genomic_position_biotypes:

				v = genomic_position_biotype[1]
				genomic_position_biotype = genomic_position_biotype[0]
				value =  v / (total_genomic_position_biotype*1.0)

				try:
					biotypes_peptides_rna_all_peptides['genomic_position_biotype'][genomic_position_biotype] += value
				except KeyError:
					biotypes_peptides_rna_all_peptides['genomic_position_biotype'][genomic_position_biotype] = value

			count_gene_level_biotypes = Counter(gene_level_biotypes_type_peptide_ribo).most_common()
			count_transcript_level_biotypes = Counter(transcript_level_biotypes_type_peptide_ribo).most_common()
			count_genomic_position_biotypes = Counter(genomic_position_biotypes_type_peptide_ribo).most_common()
			
			total_gene_level_biotype = len(gene_level_biotypes_type_peptide_ribo)
			
			for gene_level_biotype in count_gene_level_biotypes:

				v = gene_level_biotype[1]
				gene_level_biotype = gene_level_biotype[0]

				value =  v / (total_gene_level_biotype * 1.0)

				try:
					biotypes_peptides_ribo_all_peptides['gene_level_biotype'][gene_level_biotype] += value
				except KeyError:
					biotypes_peptides_ribo_all_peptides['gene_level_biotype'][gene_level_biotype] = value
				

			total_transcript_level_biotype = len(transcript_level_biotypes_type_peptide_ribo)
			
			for transcript_level_biotype in count_transcript_level_biotypes:

				v = transcript_level_biotype[1]
				transcript_level_biotype = transcript_level_biotype[0]
				value =  v / (total_transcript_level_biotype*1.0)

				try:
					biotypes_peptides_ribo_all_peptides['transcript_level_biotype'][transcript_level_biotype] += value
				except KeyError:
					biotypes_peptides_ribo_all_peptides['transcript_level_biotype'][transcript_level_biotype] = value


			total_genomic_position_biotype = len(genomic_position_biotypes_type_peptide_ribo)
			
			for genomic_position_biotype in count_genomic_position_biotypes:

				v = genomic_position_biotype[1]
				genomic_position_biotype = genomic_position_biotype[0]

				value =  v / (total_genomic_position_biotype*1.0)

				try:
					biotypes_peptides_ribo_all_peptides['genomic_position_biotype'][genomic_position_biotype] += value
				except KeyError:
					biotypes_peptides_ribo_all_peptides['genomic_position_biotype'][genomic_position_biotype] = value
				

		# https://xlsxwriter.readthedocs.io/example_pandas_multiple.html
		self.df = pd.DataFrame(data, columns=['Peptide Type', 'Peptide','Alignment', 'Strand', 'Transcript', 'gene_level_biotype', 'transcript_level_biotype', 'genomic_position_biotype', 'reads count RNA', 'reads count Ribo'])
		self.df2 = pd.DataFrame(data_resume, columns=['Peptide Type', 'Peptide','Alignment', 'Strand', 'gene_level_biotype', 'transcript_level_biotype', 'genomic_position_biotype', 'reads count RNA', 'reads count Ribo'])
		self.df3 = pd.DataFrame(data_final, columns=['Peptide Type', 'Peptide', 'gene_level_biotype', 'transcript_level_biotype', 'genomic_position_biotype', 'Total reads count RNA', 'Total reads count Ribo'])

		self.df6 = pd.DataFrame(data_unified, columns=['Peptide Type', 'Peptide','Alignment', 'Strand', 'Transcript', 'gene_level_biotype', 'transcript_level_biotype', 'genomic_position_biotype', 'ERE name', 'ERE class', 'ERE family', 'reads count RNA', 'reads count Ribo'])
		self.df7 = pd.DataFrame(data_unified_final, columns=['Peptide Type', 'Peptide', 'Consenssus Biotype', 'Total reads count RNA', 'Total reads count Ribo'])

		self.write_xls_with_all_info_biotypes()

		self.draw_biotypes(biotypes_peptides_rna, self.path_to_output_folder+'plots/biotypes/genome_annotation/biotype_by_peptide_type_reads/by_level_biotype/', True)
		self.draw_biotypes(biotypes_peptides_global, self.path_to_output_folder+'plots/biotypes/genome_annotation/global_biotypes/by_peptide_type/', False, True)
		self.draw_biotypes(biotypes_peptides_global_unified, self.path_to_output_folder+'plots/biotypes/genome_and_ERE_annotation/by_peptide_type/', False, True)
		
		print ('biotypes_peptides_all_peptides ',biotypes_peptides_all_peptides)
		self.draw_biotypes_all_peptides(biotypes_peptides_all_peptides, self.path_to_output_folder+'plots/biotypes/genome_annotation/global_biotypes/all_peptides/', False, True)
		print ('biotypes_peptides_rna_all_peptides ', biotypes_peptides_rna_all_peptides)
		self.draw_biotypes_all_peptides(biotypes_peptides_rna_all_peptides,self.path_to_output_folder+'plots/biotypes/genome_annotation/biotype_by_peptide_type_reads/all_peptides/', True)
		print ('biotypes_unified_all_peptides ', biotypes_unified_all_peptides)
		self.draw_biotypes_all_peptides(biotypes_unified_all_peptides, self.path_to_output_folder+'plots/biotypes/genome_and_ERE_annotation/all_peptides/', False)
		
		if self.mode == 'translation':
			self.draw_correlation(data)
			self.draw_biotypes(biotypes_peptides_ribo, self.path_to_output_folder+'plots/biotypes/genome_annotation/biotype_by_peptide_type_reads/by_level_biotype/', False)
			self.draw_biotypes_all_peptides(biotypes_peptides_ribo_all_peptides, self.path_to_output_folder+'plots/biotypes/genome_annotation/biotype_by_peptide_type_reads/all_peptides/', False)
			

	def draw_correlation(self, data):

		#[type_peptide, peptide, alignment, strand, 'No Annotation', gene_level_biotype, transcript_level_biotype, genomic_position_biotype, count_rna, count_ribo]
		counts_ribo = []
		counts_rna = []
		peptide_type = []
		peptides = []

		for peptide_type_, peptides_ in data.items():

			for peptide, counts in peptides_.items():
				counts_ribo.append(counts[1])
				counts_rna.append(counts[0])
				peptides.append(peptide)
			peptide_type.extend([peptide_type_]*len(peptides_))
			
		print (len(peptides))
		data = {'Type Peptide': peptide_type, 
				'Peptides': peptides,
        		'Read count Ribo': counts_ribo,
        		'Read count RNA': counts_rna} 
        		
		df_correlation = pd.DataFrame(data)
		plots.correlation(self.path_to_output_folder, self.name_exp, df_correlation)


	def draw_biotypes_2(self, biotypes_peptides, path_to_save, rna, global_ = False):

		others_non_coding = ['retained_intron', 'bidirectional_promoter_lncRNA', 'transcribed_unitary_pseudogene', 'transcribed_unprocessed_pseudogene', 'sense_overlapping','processed_pseudogene', 'unprocessed_pseudogene']
		others_protein_coding = ['IG_V_gene', 'TEC']
		others_ere = ['DNA','RC', 'RNA','Satellite','Simple_repeat','Unknown', 'Low_complexity', 'rRNA','scRNA','snRNA','srpRNA','tRNA']
		

		organisation_labels = {'Protein-coding genes':['5UTR', '3UTR', 'In_frame', 'Frameshift', 'protein_coding', 'Junctions', 'Other coding regions'], 
							'Non-coding genes':['processed_transcript', 'nonsense_mediated_decay', 'antisense', 'Exons', 'lincRNA', 'Other non-coding regions'], 
							'Protein-coding transcripts':['5UTR', '3UTR', 'In_frame', 'Frameshift', 'protein_coding', 'CDS', 'Junctions', 'Other coding regions'], 
							'Non-coding transcripts':['processed_transcript', 'nonsense_mediated_decay', 'antisense', 'Exons', 'lincRNA', 'Other non-coding regions'], 
							'Protein-coding Regions':['5UTR', '3UTR', 'In_frame', 'Frameshift', 'protein_coding', 'CDS', 'Junctions', 'Other coding regions'], 
							'Non-coding Regions':['processed_transcript', 'nonsense_mediated_decay', 'antisense', 'Exons', 'lincRNA', 'Other non-coding regions'], 
							'Intergenic Regions':['Intergenic'], 
							'Intronic Regions':['Introns'],
							'EREs':['LINE', 'LTR','Retroposon','SINE', 'Other EREs']}

		for type_peptide, level_biotypes in biotypes_peptides.items():

			for level_biotype, info_level_biotype in level_biotypes.items():

				labels_in_type_peptide = {'Protein-coding genes':{}, 'Non-coding genes': {}, 'Protein-coding Regions':{}, 'Non-coding Regions': {}, 'Protein-coding transcripts':{}, 'Non-coding transcripts': {}, 'Intergenic Regions':{}, 'Intronic Regions':{}, 'EREs':{}}
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

								if level_biotype == 'Genomic_and_ERE_Annotations_biotype':
									type_ = type_.replace('transcripts', 'Regions')
									type_ = type_.replace('genes', 'Regions')

								labels_in_type_peptide[type_][biotype] = total_biotype 
								in_ = True
								break

						if not in_:
							
							if biotype in others_non_coding:
								type_ = 'Non-coding genes'
								if level_biotype == 'transcript_level_biotype' or level_biotype == 'genomic_position_biotype':
									type_ = type_.replace('genes', 'transcripts')

								if level_biotype == 'Genomic_and_ERE_Annotations_biotype':
									type_ = type_.replace('genes', 'Regions')

								labels_in_type_peptide[type_]['Other non-coding regions'] = total_biotype 

							elif biotype in others_protein_coding:
								type_ = 'Protein-coding genes'
								if level_biotype == 'transcript_level_biotype' or level_biotype == 'genomic_position_biotype':
									type_ = type_.replace('genes', 'transcripts')

								if level_biotype == 'Genomic_and_ERE_Annotations_biotype':
									type_ = type_.replace('genes', 'Regions')

								labels_in_type_peptide[type_]['Other coding regions'] = total_biotype

							elif biotype in others_ere:
								type_ = 'EREs'
								labels_in_type_peptide[type_]['Other EREs'] = total_biotype

							else:
								if level_biotype == 'Genomic_and_ERE_Annotations_biotype':
									labels_in_type_peptide['Protein-coding Regions']['Junctions'] = total_biotype 
								else:
									labels_in_type_peptide['Protein-coding transcripts']['Junctions'] = total_biotype 

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
							plots.plot_pie(title+'rna', outer_labels, intra_labels, intra_sizes, outer_sizes, path_to_save+level_biotype+'/', self.name_exp, type_peptide+'_rna_'+level_biotype)
						else:
							plots.plot_pie(title+'ribo', outer_labels, intra_labels, intra_sizes, outer_sizes, path_to_save+level_biotype+'/', self.name_exp, type_peptide+'_ribo_'+level_biotype)
					else:
						plots.plot_pie(title, outer_labels, intra_labels, intra_sizes, outer_sizes, path_to_save, self.name_exp, type_peptide+'_'+level_biotype)
						

	def draw_biotypes_all_peptides(self, biotypes_peptides, path_to_save, rna, global_=True):

		others_non_coding = ['retained_intron', 'bidirectional_promoter_lncRNA', 'transcribed_unitary_pseudogene', 'transcribed_unprocessed_pseudogene', 'sense_overlapping','processed_pseudogene', 'unprocessed_pseudogene']
		others_protein_coding = ['IG_V_gene', 'TEC']
		others_ere = ['DNA','RC', 'RNA','Satellite','Simple_repeat','Unknown', 'Low_complexity', 'rRNA','scRNA','snRNA','srpRNA','tRNA']

		organisation_labels = {'Protein-coding genes':['5UTR', '3UTR', 'In_frame', 'Frameshift', 'protein_coding', 'Junctions', 'Other coding regions'], 
							'Non-coding genes':['processed_transcript', 'nonsense_mediated_decay', 'antisense', 'Exons', 'lincRNA', 'Other non-coding regions'], 
							'Protein-coding transcripts':['5UTR', '3UTR', 'In_frame', 'Frameshift', 'protein_coding', 'CDS', 'Junctions', 'Other coding regions'], 
							'Non-coding transcripts':['processed_transcript', 'nonsense_mediated_decay', 'antisense', 'Exons', 'lincRNA', 'Other non-coding regions'], 
							'Protein-coding Regions':['5UTR', '3UTR', 'In_frame', 'Frameshift', 'protein_coding', 'CDS', 'Junctions', 'Other coding regions'], 
							'Non-coding Regions':['processed_transcript', 'nonsense_mediated_decay', 'antisense', 'Exons', 'lincRNA', 'Other non-coding regions'], 
							'Intergenic Regions':['Intergenic'], 
							'Intronic Regions':['Introns'],
							'EREs':['LINE', 'LTR','Retroposon','SINE', 'Other EREs']}

		for level_biotype, info_level_biotype in biotypes_peptides.items():

			labels_in_type_peptide = {'Protein-coding genes':{}, 'Non-coding genes': {},'Protein-coding Regions':{}, 'Non-coding Regions': {}, 'Protein-coding transcripts':{}, 'Non-coding transcripts': {}, 'Intergenic Regions':{}, 'Intronic Regions':{}, 'EREs':{}}
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
						if level_biotype == 'Genomic_and_ERE_Annotations_biotype':
							type_ = type_.replace('transcripts', 'Regions')
							type_ = type_.replace('genes', 'Regions')

						labels_in_type_peptide[type_][biotype] = total_biotype 
						in_ = True
						break

				if not in_:
					
					if biotype in others_non_coding:
						type_ = 'Non-coding genes'
						if level_biotype == 'transcript_level_biotype' or level_biotype == 'genomic_position_biotype':
							type_ = type_.replace('genes', 'transcripts')

						if level_biotype == 'Genomic_and_ERE_Annotations_biotype':
							type_ = type_.replace('genes', 'Regions')
						labels_in_type_peptide[type_]['Other non-coding regions'] = total_biotype 

					elif biotype in others_protein_coding:
						type_ = 'Protein-coding genes'
						if level_biotype == 'transcript_level_biotype' or level_biotype == 'genomic_position_biotype':
							type_ = type_.replace('genes', 'transcripts')

						if level_biotype == 'Genomic_and_ERE_Annotations_biotype':
							type_ = type_.replace('genes', 'Regions')
						labels_in_type_peptide[type_]['Other coding regions'] = total_biotype

					elif biotype in others_ere:
						type_ = 'EREs'
						labels_in_type_peptide[type_]['Other EREs'] = total_biotype

					else:
						if level_biotype == 'Genomic_and_ERE_Annotations_biotype':
							labels_in_type_peptide['Protein-coding Regions']['Junctions'] = total_biotype 
						else:
							labels_in_type_peptide['Protein-coding transcripts']['Junctions'] = total_biotype 

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
		writer = pd.ExcelWriter(self.path_to_output_folder+'/res/Annotation_Biotypes_full.xlsx', engine='xlsxwriter')
		self.df2.to_excel(writer, sheet_name='Genomic Annotation Biotypes')
		self.df1.to_excel(writer, sheet_name='ERE Annotation Biotypes')
		self.df3.to_excel(writer, sheet_name='Genomic & ERE Annotations')
		writer.save()

	def write_xls_with_consensus_biotypes(self):
		writer = pd.ExcelWriter(self.path_to_output_folder+'/res/Annotation_Biotypes_consensus.xlsx', engine='xlsxwriter')
		self.df_total_by_position_gen.to_excel(writer, sheet_name='Genomic Annotation Biotypes')
		self.df_global_gen.to_excel(writer, sheet_name='Genomic Annotation By Region')
		self.df_final_peptide_gen.to_excel(writer, sheet_name='Genomic Annotation By peptide')
		
		self.df_total_by_position_ere.to_excel(writer, sheet_name='ERE Annotation Biotypes')
		self.df_global_ere.to_excel(writer, sheet_name='ERE Annotation By Region')
		self.df_final_peptide_ere.to_excel(writer, sheet_name='ERE Annotation By peptide')

		self.df_total_by_position_gen_ere.to_excel(writer, sheet_name='Genomic & ERE Annotations')
		self.df_consensus_annotation.to_excel(writer, sheet_name='Consensus Biotype')

		writer.save()


		