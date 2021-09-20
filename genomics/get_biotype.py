import os, logging, time, pickle, multiprocessing, _thread, csv, math, copy, pysam, copy, xlsxwriter, gc
import pandas as pd
import numpy as np
from pathos.multiprocessing import ProcessPool
import multiprocessing
import billiard as mp
import utils.useful_functions as uf
import plotting.draw_biotypes as plots
from collections import Counter, OrderedDict
from genomics.get_information_from_bed_intersection import GetInformationBEDIntersection

path_to_lib = '/'.join(os.path.abspath(__file__).split('/')[:-3])+'/lib/'

__author__ = "Maria Virginia Ruiz Cuevas"
__email__ = "maria.virginia.ruiz.cuevas@umontreal.ca"

class BiotypeAssignation:

	def __init__(self, path_to_output_folder, name_exp, mode, bam_files_list_rna, bam_files_list_ribo, order_sample_bam_files_rna, order_sample_bam_files_ribo, dev, plots):
		self.path_to_output_folder = path_to_output_folder
		self.name_exp = name_exp
		self.mode = mode
		self.plots = plots

		exists = os.path.exists(self.path_to_output_folder+'/res/BED_files/information_final_biotypes_peptides.dic') 
		exists_2 = os.path.exists(self.path_to_output_folder+'/res/BED_files/peptides_intersected_ere.dic') 
		
		if exists and exists_2:
			try:
				with open(self.path_to_output_folder+'/res/BED_files/peptides_intersected_ere.dic', 'rb') as handle:
					self.peptides_intersected_ere = pickle.load(handle)

				with open(self.path_to_output_folder+'/res/BED_files/information_final_biotypes_peptides.dic', 'rb') as handle:
					self.information_final_biotypes_peptides = pickle.load(handle)
			except ValueError:
				import pickle5
				with open(self.path_to_output_folder+'/res/BED_files/peptides_intersected_ere.dic', 'rb') as handle:
					self.peptides_intersected_ere = pickle5.load(handle)

				with open(self.path_to_output_folder+'/res/BED_files/information_final_biotypes_peptides.dic', 'rb') as handle:
					self.information_final_biotypes_peptides = pickle5.load(handle)

		else:
			self.get_info_bed_files = GetInformationBEDIntersection(path_to_output_folder)
			self.get_info_bed_files.get_information_genomic_annotation()
			self.get_info_bed_files.get_information_ERE_annotation()
			self.peptides_intersected_ere = self.get_info_bed_files.peptides_intersected_ere
			self.information_final_biotypes_peptides = self.get_info_bed_files.information_final_biotypes_peptides

		with open(path_to_lib+'ERE_info.dic', 'rb') as handle:
			self.ere_info = pickle.load(handle)

		with open(path_to_lib+'splices_information.dic', 'rb') as handle:
			self.splices_information = pickle.load(handle)

		logging.info('========== Get information from Genomic and ERE annotation : Done! ============ ')

		self.order_sample_bam_files_rna = order_sample_bam_files_rna
		self.order_sample_bam_files_ribo = order_sample_bam_files_ribo
		self.bam_files_list_rna = bam_files_list_rna
		self.bam_files_list_ribo = bam_files_list_ribo
		self.dev = dev
		
	def get_biotypes(self, info_peptide_alignments, peptides_by_type):

		logging.info('========== Getting information to define biotyping... ============ ')
		data_ere = []
		data_gen = []
		data_gen_ere = []
		self.biotype_type = set()
		self.splices_annotated = set()
		self.total_peptides = 0

		spliced = 0
		for type_peptide, peptides in peptides_by_type.items():

			for peptide in peptides:

				self.total_peptides += 1
				if True:#'AAVLEYLTAE' == peptide :#'AEKLGFAGL' == peptide:

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
							
							if '|' in position:
								chr = position.split(':')[0]
								places = position.split(':')[1].split('-')
								splices_areas = []
								in_ = False
								for place in places:
									if '|' in place:
										place_split = place.split('|')
										place_to_search = chr+':'+str(int(place_split[0])+1)+'-'+str(int(place_split[1])-1)
										try:
											in_ = place_to_search in self.splices_information[strand][chr]
											break
										except KeyError:
											pass
								if in_:
									spliced += 1
									self.splices_annotated.add(position)

							to_add_ere = []
							to_add_gen = []
							to_add_gen_ere = []

							repNames = []
							repClass_s = []
							repFamilies = []

							key_aux = peptide+'_'+position

							try:
								rep_names = self.peptides_intersected_ere[peptide][key_aux]
								
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
									#data_ere.append(to_add_ere)
									self.biotype_type.add(repClass)

							except KeyError:
								pass

							try:
								transcripts_biotypes = self.information_final_biotypes_peptides[peptide][key_aux]

								for transcript, peptide_genomic_biotypes in transcripts_biotypes.items():
									gene_level_biotype = peptide_genomic_biotypes[0]
									transcript_level_biotype = peptide_genomic_biotypes[1]
									genomic_position_biotype = peptide_genomic_biotypes[2]
										
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
									
									self.biotype_type.add(self.mod_type(genomic_position_biotype))

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
								else:
									self.biotype_type.add(genomic_position_biotype)
								
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
		columns_ere.extend(self.bam_files_list_rna)
		columns_ere.extend(self.bam_files_list_ribo)
		columns_ere.extend(['Total reads count RNA', 'Total reads count Ribo'])
		
		columns_gen = ['Peptide Type', 'Peptide','Alignment', 'MCS', 'Strand', 'Transcript', 'gene_level_biotype', 'transcript_level_biotype', 'genomic_position_biotype']
		columns_gen.extend(self.bam_files_list_rna)
		columns_gen.extend(self.bam_files_list_ribo)
		columns_gen.extend(['Total reads count RNA', 'Total reads count Ribo'])

		columns_gen_ere = ['Peptide Type', 'Peptide','Alignment', 'MCS', 'Strand', 'Transcript', 'gene_level_biotype', 'transcript_level_biotype', 'genomic_position_biotype', 'ERE name', 'ERE class', 'ERE family']
		columns_gen_ere.extend(self.bam_files_list_rna)
		columns_gen_ere.extend(self.bam_files_list_ribo)
		columns_gen_ere.extend(['Total reads count RNA', 'Total reads count Ribo'])

		#self.df1 = pd.DataFrame(data_ere, columns = columns_ere)
		self.df2 = pd.DataFrame(data_gen, columns = columns_gen)
		self.df3 = pd.DataFrame(data_gen_ere, columns = columns_gen_ere)
		
		data_gen_ere = []
		data_gen = []
		data_ere = []

		logging.info('========== Getting information to define biotyping : Done! ============ ')

	def mod_type(self, type_):
		if '-' in type_:
			if 'Non_coding' in type_ :
				type_ = 'Non_coding Junctions'
			else:
				type_ = 'Junctions'
		return type_

	def get_global_annotation(self):

		# ERE Full
		groupby_columns = ['Peptide Type', 'Peptide','Alignment', 'Strand', 'ERE name', 'ERE class', 'ERE family']
		bam_files_columns = self.bam_files_list_rna
		bam_files_columns.extend(self.bam_files_list_ribo)
		bam_files_columns.extend(['Total reads count RNA', 'Total reads count Ribo'])
		
		#self.df_total_by_position_ere = self.df1.groupby(groupby_columns)[bam_files_columns].sum().reset_index()
		#self.df_total_by_position_ere = pd.DataFrame(self.df_total_by_position_ere, columns = groupby_columns+bam_files_columns)
		
		# Genomic Annotation Full
		groupby_columns = ['Peptide Type', 'Peptide','Alignment', 'Strand', 'Transcript', 'gene_level_biotype', 'transcript_level_biotype', 'genomic_position_biotype']
		self.df_total_by_position_gen = self.df2.groupby(groupby_columns)[bam_files_columns].sum().reset_index()
		self.df_total_by_position_gen = pd.DataFrame(self.df_total_by_position_gen, columns = groupby_columns+bam_files_columns)
		del self.df2

		# Genomic and ERE Annotation Full
		groupby_columns = ['Peptide Type', 'Peptide','Alignment', 'Strand', 'Transcript', 'gene_level_biotype', 'transcript_level_biotype', 'genomic_position_biotype', 'ERE name', 'ERE class', 'ERE family']
		self.df_total_by_position_gen_ere = self.df3.groupby(groupby_columns)[bam_files_columns].sum().reset_index()
		self.df_total_by_position_gen_ere = pd.DataFrame(self.df_total_by_position_gen_ere, columns = groupby_columns+bam_files_columns)
		
		shape_df3 = self.df3.shape
		cols = shape_df3[1]

		if self.total_peptides >= 2000 and cols >= 2000 :
			path = self.path_to_output_folder+'/res/full_info_biotypes/1_Genomic_and_ERE_Annotations_Full.csv'
			self.df3.to_csv(path, index=False)
			del self.df3

		self.get_global_annotation_gen(bam_files_columns)
		self.get_genomic_and_ere_annotation(bam_files_columns)

		logging.info('========== Global genomic and ERE annotation : Done! ============ ')
		
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
		data_global_annotation_gen = []
		df_global_gen_construction = []
		#with open(self.path_to_output_folder+'alignments/df_global_gen.dic', 'wb') as handle:
		#	pickle.dump(self.df_global_gen, handle, protocol=pickle.HIGHEST_PROTOCOL)

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

	def get_information_mixed_genomic_ERE(self, df_final_gen_sum, group_by_gen_ere):

		# This method gets the information for each alignment and gives the 
		# weighted percentage accordingly each biotype for a position. 
		# ERE and Genomic are mixed.

		peptides_visited = {}
		peptides_visited_alignment = {}
		only_annotated_splices = {}

		def insert_dic(biotype, sum_percentage, in_):

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
			
			if in_:
				# This is only to keep a trace of those alignments that are annotated splices. 
				try:
					info_to_subtract = only_annotated_splices[(peptide_type, peptide)]
					try:
						info_to_subtract[biotype] += sum_percentage
					except KeyError:
						info_to_subtract[biotype] = sum_percentage
				except KeyError:
					dic2 = {}
					dic2[biotype] = sum_percentage
					only_annotated_splices[(peptide_type, peptide)] = dic2
			try:
				biotypes = peptides_visited_alignment[key][0]
				try:
					biotypes[biotype] += sum_percentage
				except KeyError:
					biotypes[biotype] = sum_percentage
			except KeyError:
				dic = {}
				dic[biotype] = sum_percentage
				peptides_visited_alignment[key] = [dic]			

		for key, items in group_by_gen_ere:
			peptide_type = key[0]
			peptide = key[1]
			alignment = key[2]

			in_ = True

			if '|' in alignment:
				# Gets the information if the splice alignment is annotated or not, in order to filter and keep those that 
				# are annotated.
				in_ = alignment in self.splices_annotated
			
			count_bamfiles = df_final_gen_sum.get_group((peptide_type, peptide)).values.tolist()[0]
			count_bamfiles_alignment = group_by_gen_ere.get_group((key)).values.tolist()[0][11:]
			
			type_gen_aux = list(filter(lambda a:a != '', items['genomic_position_biotype'].values.tolist()))
			type_ere_aux = list(filter(lambda a:a != '', items['ERE class'].values.tolist()))

			type_ere = []
			type_gen = []

			if 'In_frame' not in type_gen_aux:
				for type_ in type_ere_aux:
					types = type_.split(',')
					type_ere.extend(types)

				for type_ in type_gen_aux:
					types = type_.split(',')
					type_gen.extend(types)
			else:
				type_gen = ['In_frame']

			type_ere = list(type_ere)
			type_gen = list(type_gen)

			total_items = len(type_gen) + len(type_ere)
			sum_percentage = 1 / (total_items*1.0)
			
			for genomic_position_biotype in type_gen:
				if genomic_position_biotype != '' :
					insert_dic(genomic_position_biotype, sum_percentage, in_)
			
			for ere_class in type_ere:
				if ere_class != '' :
					insert_dic(ere_class, sum_percentage, in_)
			
			if len(peptides_visited[(peptide_type, peptide)]) == 1:
				peptides_visited[(peptide_type, peptide)].append(count_bamfiles)

			if len(peptides_visited_alignment[key]) == 1:
				peptides_visited_alignment[key].append(count_bamfiles_alignment)

		return peptides_visited_alignment, peptides_visited, only_annotated_splices

	def genome_annotations(self, general_info_peptides, only_annotated_splices, bam_files_columns):

		self.counts_reads_rna_ribo = {} # This is for plot correlation if Ribosome Profiling information is available
		df_consensus_annotation = [] # General Gen & ERE Biotype
		peptides_absent_sample_group = {}
		biotype_info_all_alignments_annotated = {} # Taking into account all the alignments
		biotype_info_only_alignments_annotated = {} # Not taking into account the splices sites non annotated
		biotypes_by_peptide_genome_explained = []
		biotypes_by_peptide_type = {}
		biotypes_all_peptides = {}
		biotype_type = list(self.biotype_type)

		for key, information_peptide in general_info_peptides.items():
			
			to_add = []
			to_add_aux = []
			peptide_type = key[0]
			peptide = key[1]
			
			total_alignments_pep = sum(general_info_peptides[(peptide_type, peptide)][0].values())
			biotypes_all_alignments = general_info_peptides[(peptide_type, peptide)][0]
			count_bamfiles = general_info_peptides[(peptide_type, peptide)][1]

			try:
				biotypes_annotated_alignments = only_annotated_splices[(peptide_type, peptide)]
				total_alignments_spliced_annotated = sum(only_annotated_splices[(peptide_type, peptide)].values())
			except KeyError:
				biotypes_annotated_alignments = ''
				total_alignments_spliced_annotated = 0
			
			count_rna = count_bamfiles[-2]
			count_ribo = count_bamfiles[-1]

			count_string = []
			best_guess = ''
			count_bamfiles_ = count_bamfiles[:-2]
			to_add_aux.extend(key)

			bios = [0]*(len(biotype_type)+1)
			
			if 0 in count_bamfiles_:
				indices = [i for i, x in enumerate(count_bamfiles_) if x == 0]
				peptides_absent_sample_group[peptide] = indices

			if biotypes_annotated_alignments != '':
				in_frame = False
				for type_ in sorted(biotypes_annotated_alignments, key=biotypes_annotated_alignments.get, reverse=True):
					# if 'In_frame' in biotypes_all_alignments.keys():
					# 	in_frame = True
					# 	count = 1
					# 	type_ = self.mod_type(type_)
					# else:
					count = biotypes_annotated_alignments[type_]
					count = count/(total_alignments_spliced_annotated*1.0)
					type_ = self.mod_type(type_)

					try:
						type_in_dic  = biotype_info_only_alignments_annotated['Genome']
						try:
							types_dic = type_in_dic[peptide_type]
							try:
								types_dic[type_]+= count
							except KeyError:
								types_dic[type_] = count
						except KeyError:
							type_in_dic[peptide_type]= {type_ : count}
					except KeyError:
						peptide_type_dic = {}
						peptide_type_dic[peptide_type] = {type_ : count}
						biotype_info_only_alignments_annotated['Genome'] = peptide_type_dic

					if in_frame:
						break

			for type_ in sorted(biotypes_all_alignments, key=biotypes_all_alignments.get, reverse=True):

				if type_ == 'In_frame':
					best_guess = 'In_frame'
				
				count = biotypes_all_alignments[type_]
				count = count/(total_alignments_pep*1.0)
				
				percentage = round(count*100,2)
				count_string.append(type_+': '+str(percentage)+'%')

				# Plots
				if self.plots:
					peptide_type_aux = peptide_type.split(',')

					for pep_type in peptide_type_aux:
						
						try:
							types = biotypes_by_peptide_type[pep_type]
							try:
								types[type_] += count
							except KeyError:
								types[type_] = count
						except KeyError:
							biotypes_by_peptide_type[pep_type] = {type_ : count}

					try:
						biotypes_all_peptides[type_] += count
					except KeyError:
						biotypes_all_peptides[type_] = count


				type_ = self.mod_type(type_)
				type_index = biotype_type.index(type_)
				bios[type_index] += percentage
				
				try:
					type_in_dic  = biotype_info_all_alignments_annotated['Genome']
					try:
						type_in_dic[type_]+= count
					except KeyError:
						type_in_dic[type_] = count
				except KeyError:
					biotype_info_all_alignments_annotated['Genome'] = {type_ : count}

			try:
				self.counts_reads_rna_ribo[peptide_type][peptide] = [count_rna, count_ribo]
			except KeyError:
				dic = {}
				dic[peptide] = [count_rna, count_ribo]
				self.counts_reads_rna_ribo[peptide_type] = dic
			
			consensus = ' - '.join(count_string)
			
			if best_guess == '':
				best_guess = consensus
				max_best_guess = max(bios)
				indexes = [i for i,x in enumerate(bios) if x==max_best_guess]
				guessess = map(biotype_type.__getitem__, indexes)
				best_guess = ','.join(guessess)
				bios[-1] = best_guess

			to_add_aux.extend(bios)
			to_add = [peptide_type, peptide, consensus, best_guess]
			to_add.extend(count_bamfiles)
			df_consensus_annotation.append(to_add)
			if self.dev:
				biotypes_by_peptide_genome_explained.append(to_add_aux)

		groupby_columns = ['Peptide Type', 'Peptide', 'Biotype Frequencies', 'Best Guess']
		self.df_consensus_annotation = pd.DataFrame(df_consensus_annotation, columns = groupby_columns+bam_files_columns)

		shape_df_consensus_annotation = self.df_consensus_annotation.shape
		cols = shape_df_consensus_annotation[1]

		if self.total_peptides >= 2000 and cols >= 2000 :
			path = self.path_to_output_folder+'/res/summary_info_biotypes/1_General_Gen_and_ERE_Biotype_Consensus.csv'
			self.df_consensus_annotation.to_csv(path, index=False)
			del self.df_consensus_annotation

		df_consensus_annotation = []
		
		if self.plots:
			logging.info('========== Plots ============ ')
			
			plots.draw_biotypes(biotypes_by_peptide_type, self.path_to_output_folder+'plots/biotypes/genome_and_ERE_annotation/by_peptide_type/', False, False, self.name_exp)
			logging.info('========== biotypes_by_peptide_type : Genome & ERE annotations : Done! ============ ')
			
			plots.draw_biotypes(biotypes_all_peptides, self.path_to_output_folder+'plots/biotypes/genome_and_ERE_annotation/all_peptides/', True, False, self.name_exp)
			logging.info('========== biotypes_all_peptides : Genome & ERE annotations : Done! ============ ')
			
			logging.info('========== Plots : Done! ============ ')
			
		if self.mode == 'translation':
			plots.draw_correlation(counts_reads_rna_ribo, self.name_exp, self.path_to_output_folder)

		
		if self.dev:
			groupby_columns = ['Peptide Type', 'Peptide']
			cols_names = groupby_columns+biotype_type
			cols_names.append('Best Guess')
			self.biotypes_by_peptide_genome_explained = pd.DataFrame(biotypes_by_peptide_genome_explained, columns = cols_names)

			path = self.path_to_output_folder+'/res/full_info_biotypes/biotypes_by_peptide_genome_explained.csv'
			self.biotypes_by_peptide_genome_explained.to_csv(path, index=False)

			del self.biotypes_by_peptide_genome_explained
			
			with open(self.path_to_output_folder+'alignments/biotypes_by_peptide_genome_explained.list', 'wb') as handle:
				pickle.dump(biotypes_by_peptide_genome_explained, handle, protocol=pickle.HIGHEST_PROTOCOL)

			biotypes_by_peptide_genome_explained = []
			gc.collect()

		return peptides_absent_sample_group, biotype_info_all_alignments_annotated, biotype_info_only_alignments_annotated 

	def transcription_annotations(self, info_peptides_by_alignment, general_info_peptides, indexes_group, indexes_by_group, split_bams_files, biotype_info_all_alignments_annotated, biotype_info_only_alignments_annotated, bam_files_columns):
		
		info_peptides_alignments_by_sample = {}
		info_peptides_alignments_by_sample_group = {}
		biotypes_all_peptides_group_samples = {}
		biotypes_by_peptide_type_group_samples = {}
		biotypes_all_peptides_type_group_samples_all = {}
		df_consensus_annotation_full = []

		keys = list(info_peptides_by_alignment.keys())
		values = list(info_peptides_by_alignment.values())
		
		def fill_information(res):
			to_add = res[0]
			info_peptides_alignments_by_sample_aux = res[1]
			
			# index, type_, value, key_group, value_group
			for index, element in enumerate(info_peptides_alignments_by_sample_aux):
				if index == 0:
					key = element
					peptide_type = key[0]
					peptide_type_aux = peptide_type.split(',')
				elif index == 1:
					count_bamfiles = element
				else:
					index = element[0]
					type_ = element[1]
					value = element[2]
					key_group = element[3]
					value_group = element[4]
					
					try:
						info_peptide = info_peptides_alignments_by_sample[key]
						try:
							info_peptide[index][type_] += value
						except TypeError:
							dic_aux = {type_ : value}
							info_peptides_alignments_by_sample[key][index] = dic_aux
						except KeyError:
							info_peptide[index][type_] = value
							
					except KeyError:
						dic_aux = [0]*len(count_bamfiles)
						dic_aux[index] = {type_ : value}
						info_peptides_alignments_by_sample[key] = dic_aux

					try:
						info_peptide = info_peptides_alignments_by_sample_group[key]
						try:
							info_peptide[key_group][type_] += value_group
						except TypeError:
							dic_aux = {type_ : value_group}
							info_peptides_alignments_by_sample_group[key][key_group] = dic_aux
						except KeyError:
							info_peptide[key_group][type_] = value_group
					except KeyError:
						dic_aux = copy.deepcopy(split_bams_files)
						dic_aux[key_group] = {type_ : value_group}
						info_peptides_alignments_by_sample_group[key] = dic_aux

					if key_group != 'Total reads count RNA' and key_group != 'Total reads count Ribo':
						
						if self.plots:
							try:
								types = biotypes_all_peptides_group_samples[key_group]
								try:
									types[type_] += value_group

								except KeyError:
									biotypes_all_peptides_group_samples[key_group][type_] = value_group
									
							except KeyError:
								biotypes_all_peptides_group_samples[key_group] = {type_ : value_group}

							for pep_type in peptide_type_aux:
								try:
									peptide_types = biotypes_by_peptide_type_group_samples[key_group]
									try:
										types = peptide_types[pep_type]
										try:
											types[type_] += value_group
										except KeyError:
											types[type_] = value_group
									except KeyError:
										peptide_types[pep_type] = {type_ : value_group}
								except KeyError:
									dic_aux = {}
									dic_aux[pep_type] = {type_ : value_group}
									biotypes_by_peptide_type_group_samples[key_group] = dic_aux

					if key_group == 'Total reads count RNA':

						if self.plots:
							try:
								biotypes_all_peptides_type_group_samples_all[type_] += value_group

							except KeyError:
								biotypes_all_peptides_type_group_samples_all[type_] = value_group

						type_ = self.mod_type(type_)

						if value != 0:
							try:
								type_in_dic  = biotype_info_all_alignments_annotated['All']
								try:
									type_in_dic[type_]+= value
								except KeyError:
									type_in_dic[type_] = value
							except KeyError:
								biotype_info_all_alignments_annotated['All'] = {type_ : value}

							try:
								type_in_dic  = biotype_info_only_alignments_annotated['All']
								try:
									types_dic = type_in_dic[peptide_type]
									try:
										types_dic[type_]+= value
									except KeyError:
										types_dic[type_] = value
								except KeyError:
									type_in_dic[peptide_type]= {type_ : value}
							except KeyError:
								peptide_type_dic = {}
								peptide_type_dic[peptide_type] = {type_ : value}
								biotype_info_only_alignments_annotated['All'] = peptide_type_dic

			df_consensus_annotation_full.append(to_add)
		
		indexes_group_to_send = [indexes_group]*len(keys)
		indexes_by_group_to_send = [indexes_by_group]*len(keys)
			
		if len(keys) > 50000:
			# cont = 0
			# nodes = 5
			# for i in range(0,len(keys),nodes):
			# 	#pool = mp.Pool(processes = 1, maxtasksperchild=1000)
			# 	print (i)
			# 	pool = ProcessPool(nodes=nodes)
			# 	cont += nodes
			# 	#chunk = tasks[i:cont]
			# 	#results = pool.starmap(self.get_information_transcription, chunk)
			# 	results = pool.map(self.get_information_transcription, keys[i:cont], values[i:cont], indexes_group_to_send[i:cont], indexes_by_group_to_send[i:cont])
			# 	for res in results:
			# 		fill_information(res)

			# 	print (i)
			# 	pool.close()
			# 	pool.join()
			# 	pool.clear()
			for i in range(0,len(keys)):
				results = self.get_information_transcription (keys[i], values[i], indexes_group, indexes_by_group)
				fill_information(results)				
		else:
			
			#tasks = [*zip(keys, values, indexes_group_to_send, indexes_by_group_to_send)]
			nodes = multiprocessing.cpu_count()
			#pool = mp.Pool(nodes)
			pool = ProcessPool(nodes=nodes)
			#results = pool.starmap(self.get_information_transcription, tasks)
			results = pool.map(self.get_information_transcription, keys, values, indexes_group_to_send, indexes_by_group_to_send)

			for res in results:
				fill_information(res)

			pool.close()
			pool.join()
			pool.clear()

		groupby_columns = ['Peptide Type', 'Peptide', 'Alignment', 'Strand', 'Consenssus']
		self.df_consensus_annotation_full = pd.DataFrame(df_consensus_annotation_full, columns = groupby_columns+bam_files_columns)
		
		shape_df_consensus_annotation_full = self.df_consensus_annotation_full.shape
		cols = shape_df_consensus_annotation_full[1]

		try:
			l = len(self.df3)
		except AttributeError:
			path = self.path_to_output_folder+'/res/full_info_biotypes/3_Genomic_and_ERE_Anno_by_Region_Full.csv'
			self.df_consensus_annotation_full.to_csv(path, index=False)
			del self.df_consensus_annotation_full 

		df_consensus_annotation_full = []
		info_peptides_by_alignment = []
		gc.collect()

		if self.plots:
			logging.info('========== Plots ============ ')

			plots.draw_biotypes(biotypes_all_peptides_group_samples, self.path_to_output_folder+'plots/biotypes/biotype_by_sample_group/all_peptides/', True, True, self.name_exp)
			logging.info('========== biotypes_all_peptides_group_samples : Transcription levels : Done! ============ ')
			
			plots.draw_biotypes(biotypes_by_peptide_type_group_samples, self.path_to_output_folder+'plots/biotypes/biotype_by_sample_group/by_peptide_type/', False, True, self.name_exp)
			logging.info('========== biotypes_by_peptide_type_group_samples : Transcription levels : Done! ============ ')
			
			plots.draw_biotypes(biotypes_all_peptides_type_group_samples_all, self.path_to_output_folder+'plots/biotypes/biotype_by_sample_group/all_peptides/', True, False, self.name_exp)
			logging.info('========== biotypes_all_peptides_type_group_samples_all : Transcription levels : Done! ============ ')
			
			logging.info('========== Plots : Done! ============ ')

		return info_peptides_alignments_by_sample, info_peptides_alignments_by_sample_group

	def get_information_transcription(self, key, information_peptide, indexes_group, indexes_by_group):
		
		to_add = []
		peptide_type = key[0]
		peptide = key[1]
		alignment = key[2]
		strand = key[3]

		total_alignments_pep = sum(information_peptide[0].values())
		biotypes_by_alignment = information_peptide[0]
		
		count_bamfiles = information_peptide[1]
		count_bamfiles_total = self.general_info_peptides[(peptide_type, peptide)][1]
		
		count_rna = count_bamfiles[-2]
		count_ribo = count_bamfiles[-1]

		count_string = []
		to_add = [peptide_type, peptide, alignment, strand]
		
		info_peptides_alignments_by_sample = [(peptide_type, peptide), count_bamfiles]
		
		for type_ in sorted(biotypes_by_alignment, key=biotypes_by_alignment.get, reverse=True):
			count = biotypes_by_alignment[type_]
			count = count/(total_alignments_pep*1.0)
			
			percentage = str(round(count*100,2))
			count_string.append(type_+': '+percentage+'%')
			
			for index, bamfile in enumerate(count_bamfiles):
				value_in_bam_file = bamfile
				value_total_bam_file = count_bamfiles_total[index]
				key_group = indexes_group[index]
				total_group_alignment = np.array(count_bamfiles_total)[indexes_by_group[key_group]].sum()
				try:
					value = count*(value_in_bam_file/value_total_bam_file)
				except ZeroDivisionError:
					value = math.log(1,10)

				try:
					value_group = value*(value_total_bam_file/total_group_alignment)
				except ZeroDivisionError:
					value_group = 0
				if value == 0:
					value_group = 0

				info_peptides_alignments_by_sample.append((index, type_, value, key_group, value_group))

		consensus = ' - '.join(count_string)
		to_add.append(consensus)
		to_add.extend(count_bamfiles)

		return to_add, info_peptides_alignments_by_sample

	def process_information_by_peptides_in_samples(self, info_peptides_alignments_by_sample, biotype_info_all_alignments_annotated, biotype_info_only_alignments_annotated, bam_files_columns):

		df_consensus_annotation_final = []
		biotypes_by_peptide_sample_explained = []
		biotype_type = list(self.biotype_type)

		for key, information_peptide in info_peptides_alignments_by_sample.items():
			to_add_aux = []
			to_add_aux.extend(key)
			peptide_type = key[0]

			for index, bam_file_dic in enumerate(information_peptide):
				count_string = []
				to_add = []
				sample = bam_files_columns[index]
				to_add.append(sample)
				to_add.extend(key)
				bios = [0]*(len(biotype_type)+1)

				for type_ in sorted(bam_file_dic, key=bam_file_dic.get, reverse=True):
					count = bam_file_dic[type_]
					type_ = self.mod_type(type_)
					type_index = biotype_type.index(type_)
					
					if count > 0:
						percentage = round(count*100,2)
						count_string.append(type_+': '+str(percentage)+'%')
						bios[type_index] += percentage

						if sample != 'Total reads count RNA':
							try:
								type_in_dic  = biotype_info_all_alignments_annotated[sample]
								try:
									type_in_dic[type_] += count
								except KeyError:
									type_in_dic[type_] = count
							except KeyError:
								biotype_info_all_alignments_annotated[sample] = {type_ : count}

							try:
								type_in_dic  = biotype_info_only_alignments_annotated[sample]
								try:
									types_dic = type_in_dic[peptide_type]
									try:
										types_dic[type_]+= count
									except KeyError:
										types_dic[type_] = count
								except KeyError:
									type_in_dic[peptide_type]= {type_ : count}
							except KeyError:
								peptide_type_dic = {}
								peptide_type_dic[peptide_type] = {type_ : count}
								biotype_info_only_alignments_annotated[sample] = peptide_type_dic

				if sum(bios) != 0:
					max_best_guess = max(bios)
					indexes = [i for i,x in enumerate(bios) if x==max_best_guess]
					guessess = map(biotype_type.__getitem__, indexes)
					best_guess = ','.join(guessess)
					bios[-1] = best_guess
				else:
					bios[-1] = 'NA'
				to_add.extend(bios)
				consensus = ' - '.join(count_string)
				to_add_aux.append(consensus)
				if self.dev:
					biotypes_by_peptide_sample_explained.append(to_add)
			df_consensus_annotation_final.append(to_add_aux)

		with open(self.path_to_output_folder+'alignments/biotype_info_all_alignments_annotated.dic', 'wb') as handle:
			pickle.dump(biotype_info_all_alignments_annotated, handle, protocol=pickle.HIGHEST_PROTOCOL)
		biotype_info_all_alignments_annotated = []
		gc.collect()

		with open(self.path_to_output_folder+'alignments/biotype_info_only_alignments_annotated.dic', 'wb') as handle:
			pickle.dump(biotype_info_only_alignments_annotated, handle, protocol=pickle.HIGHEST_PROTOCOL)
		biotype_info_only_alignments_annotated = []
		gc.collect()

		groupby_columns = ['Peptide Type', 'Peptide']
		self.df_consensus_annotation_full_final = pd.DataFrame(df_consensus_annotation_final, columns = groupby_columns+bam_files_columns)
		
		try:
			l = len(self.df_consensus_annotation)
		except AttributeError:
			path = self.path_to_output_folder+'/res/summary_info_biotypes/2_Sample_Gen_and_ERE_Biotype_Consensus.csv'
			self.df_consensus_annotation_full_final.to_csv(path, index=False)
			del self.df_consensus_annotation_full_final

		df_consensus_annotation_final = []

		if self.dev:
			groupby_columns = ['Sample', 'Peptide Type', 'Peptide']
			cols_names = groupby_columns+biotype_type
			cols_names.append('Best Guess')
			self.biotypes_by_peptide_sample_explained = pd.DataFrame(biotypes_by_peptide_sample_explained, columns = cols_names)
			path = self.path_to_output_folder+'/res/full_info_biotypes/biotypes_by_peptide_sample_explained.csv'
			self.biotypes_by_peptide_sample_explained.to_csv(path, index=False)
			del self.biotypes_by_peptide_sample_explained
			biotypes_by_peptide_sample_explained = []
			gc.collect()

	def get_genomic_and_ere_annotation(self, bam_files_columns):

		logging.info('========== Init: Genomic and ERE annotation ============ ')
		
		groupby_columns = ['Peptide Type', 'Peptide', 'Alignment', 'Strand']
		group_by_gen_ere = self.df_total_by_position_gen_ere.groupby(groupby_columns)

		shape_df_total_by_position_gen_ere = self.df_total_by_position_gen_ere.shape
		cols = shape_df_total_by_position_gen_ere[1]

		try:
			l = len(self.df3)
		except AttributeError:
			path = self.path_to_output_folder+'/res/full_info_biotypes/2_Genomic_and_ERE_Annotations_Summary_Full.csv'
			self.df_total_by_position_gen_ere.to_csv(path, index=False)
			del self.df_total_by_position_gen_ere 

		groupby_columns = ['Peptide Type', 'Peptide']
		df_final_gen_sum = self.df_global_gen.groupby(groupby_columns)[bam_files_columns].sum()
		df_final_gen_sum = df_final_gen_sum.groupby(groupby_columns)
		
		del self.df_global_gen
		
		# Gets the information for each peptide with the information for each alignment,
		# and for all the alingments.

		res = self.get_information_mixed_genomic_ERE(df_final_gen_sum, group_by_gen_ere)
		info_peptides_by_alignment = res[0] 
		self.general_info_peptides = res[1] 
		only_annotated_splices = res[2]
		logging.info('========== Get information from genomic and ere intersection : Done!  ============ ')

		res = self.genome_annotations(self.general_info_peptides, only_annotated_splices, bam_files_columns)
		peptides_absent_sample_group = res[0] 
		biotype_info_all_alignments_annotated = res[1] 
		biotype_info_only_alignments_annotated = res[2] 
		
		logging.info('========== Get Biotype based on absence of transcription (genome based) : Done!  ============ ')

		res = self.get_info_group_samples(bam_files_columns)
		split_bams_files = res[0] 
		indexes_group = res[1] 
		order_columns = res[2] 
		indexes_by_group = res[3]

		res = self.transcription_annotations(info_peptides_by_alignment, self.general_info_peptides, indexes_group, indexes_by_group, split_bams_files, biotype_info_all_alignments_annotated, biotype_info_only_alignments_annotated, bam_files_columns)
		info_peptides_alignments_by_sample = res[0] 
		info_peptides_alignments_by_sample_group = res[1]

		logging.info('========== Get Biotype based on transcription : Done!  ============ ')

		self.get_biotype_by_group_sample(order_columns, info_peptides_alignments_by_sample_group)
		
		logging.info('========== Writting full biotyping in xls files ============ ')

		self.write_xls_with_all_info_biotypes()
		gc.collect()

		logging.info('========== Writting full biotyping in xls files : Done! ============ ')

		self.process_information_by_peptides_in_samples(info_peptides_alignments_by_sample, biotype_info_all_alignments_annotated, biotype_info_only_alignments_annotated, bam_files_columns)
	
		logging.info('========== Writting consensus biotyping in xls files ============ ')

		self.write_xls_with_consensus_biotypes()
		gc.collect()

		logging.info('========== Writting consensus biotyping in xls files : Done! ============ ')

		logging.info('========== Fini: Genomic and ERE annotation : Done! ============ ')

	def get_info_group_samples(self, bam_files_columns):

		split_bams_files = {} 
		indexes_group = {}
		order_columns = []
		indexes_by_group = {}

		for group, samples in self.order_sample_bam_files_rna.items():
			order_columns.append(group)
			for sample in samples:
				index = bam_files_columns.index(sample)
				indexes_group[index] = group
				split_bams_files[group] = {}
				try:
					indexes_by_group[group].append(index)
				except KeyError:
					indexes_by_group[group] = [index]

		for group, samples in self.order_sample_bam_files_ribo.items():
			order_columns.append(group)
			for sample in samples:
				index = bam_files_columns.index(sample)
				indexes_group[index] = group
				split_bams_files[group] = {}
				try:
					indexes_by_group[group].append(index)
				except KeyError:
					indexes_by_group[group] = [index]

		indexes_group[len(bam_files_columns)-2] = 'Total reads count RNA'
		indexes_group[len(bam_files_columns)-1] = 'Total reads count Ribo'
		indexes_by_group['Total reads count RNA'] = [len(bam_files_columns)-2]
		indexes_by_group['Total reads count Ribo'] = [len(bam_files_columns)-1]

		split_bams_files['Total reads count RNA'] = {} 
		split_bams_files['Total reads count Ribo'] = {}
		order_columns.append('Total reads count RNA')
		order_columns.append('Total reads count Ribo')

		return split_bams_files, indexes_group, order_columns, indexes_by_group

	def get_biotype_by_group_sample(self, order_columns, peptides_visited_sample_group):

		df_consensus_annotation_final_sample = []
		
		for key, information_peptide in peptides_visited_sample_group.items():
			to_add_aux = []
			to_add_aux.extend(key)
			for group in order_columns:
				group_dic = information_peptide[group]
				count_string = []
				for type_ in sorted(group_dic, key=group_dic.get, reverse=True):
					count = group_dic[type_]
					if count > 0:
						percentage = str(round(count*100,2))
						count_string.append(type_+': '+percentage+'%')
				consensus = ' - '.join(count_string)
				to_add_aux.append(consensus)

			df_consensus_annotation_final_sample.append(to_add_aux)
		groupby_columns = ['Peptide Type', 'Peptide']
		self.df_consensus_annotation_final_sample = pd.DataFrame(df_consensus_annotation_final_sample, columns = groupby_columns+order_columns)
		df_consensus_annotation_final_sample = []

		try:
			l = len(self.df_consensus_annotation)
		except AttributeError:
			path = self.path_to_output_folder+'/res/summary_info_biotypes/3_Group_Samples_Gen_and_ERE_Biotype_Consensus.csv'
			self.df_consensus_annotation_final_sample.to_csv(path, index=False)
			del df_consensus_annotation_final_sample


	def write_xls_with_all_info_biotypes(self):

		exist_df3 = True
		try:
			l = len(self.df3)
		except AttributeError:
			exist_df3 = False

		if  exist_df3 :

			writer = pd.ExcelWriter(self.path_to_output_folder+'/res/full_info_biotypes/Annotation_Biotypes_full_info.xlsx', engine='xlsxwriter')
			writer.book.use_zip64()
			
			self.df3.to_excel(writer, sheet_name='Genomic & ERE Annotations')
			worksheet1 = writer.sheets['Genomic & ERE Annotations']
			worksheet1.set_tab_color('purple')
			del self.df3

			
			self.df_total_by_position_gen_ere.to_excel(writer, sheet_name='Genomic & ERE Annotations_')
			worksheet1 = writer.sheets['Genomic & ERE Annotations_']
			worksheet1.set_tab_color('purple')
			del self.df_total_by_position_gen_ere 

			
			self.df_consensus_annotation_full.to_excel(writer, sheet_name='Genomic & ERE Anno. By Region')
			worksheet1 = writer.sheets['Genomic & ERE Anno. By Region']
			worksheet1.set_tab_color('purple')
			del self.df_consensus_annotation_full 

			#https://xlsxwriter.readthedocs.io/working_with_colors.html
			writer.save()

	def write_xls_info_biotypes_explained(self):
		
		if (len(self.biotypes_by_peptide_sample_explained) or len(self.biotypes_by_peptide_sample_explained) )  < 1048576:
			writer = pd.ExcelWriter(self.path_to_output_folder+'/res/full_info_biotypes/biotypes_by_peptide_sample_explained.xlsx', engine='xlsxwriter')
			writer.book.use_zip64()
			self.biotypes_by_peptide_sample_explained.to_excel(writer, sheet_name='Biotypes Sample Explained')
			self.biotypes_by_peptide_genome_explained.to_excel(writer, sheet_name='Biotypes Genome Explained')
			
			writer.save()
		else:
			logging.info('Failing to write xls files: information overpasses the limit for rows into a xls file. ')
			logging.info('========== Writting individual information of biotyping in csv files ============ ')
			path = self.path_to_output_folder+'/res/full_info_biotypes/biotypes_by_peptide_sample_explained.csv'
			self.biotypes_by_peptide_sample_explained.to_csv(path, index=False)
			path = self.path_to_output_folder+'/res/full_info_biotypes/biotypes_by_peptide_genome_explained.csv'
			self.biotypes_by_peptide_genome_explained.to_csv(path, index=False)

	def write_xls_with_consensus_biotypes(self):

		exist_df_consensus_annotation = True
		try:
			l = len(self.df_consensus_annotation)
		except AttributeError:
			exist_df_consensus_annotation = False

		if  exist_df_consensus_annotation:

			writer = pd.ExcelWriter(self.path_to_output_folder+'/res/summary_info_biotypes/Annotation_Biotypes_consensus.xlsx', engine='xlsxwriter')
			writer.book.use_zip64()

			
			self.df_consensus_annotation.to_excel(writer, sheet_name='General Gen & ERE Biotype')
			worksheet1 = writer.sheets['General Gen & ERE Biotype']
			worksheet1.set_tab_color('pink')
			del self.df_consensus_annotation

		
			self.df_consensus_annotation_full_final.to_excel(writer, sheet_name='Sample Gen & ERE Biotype')
			worksheet1 = writer.sheets['Sample Gen & ERE Biotype']
			worksheet1.set_tab_color('purple')
			del self.df_consensus_annotation_full_final

		
			self.df_consensus_annotation_final_sample.to_excel(writer, sheet_name='Group Samples Gen & ERE Biotype')
			worksheet1 = writer.sheets['Group Samples Gen & ERE Biotype']
			worksheet1.set_tab_color('navy')
			del self.df_consensus_annotation_final_sample

			writer.save()
			
			





