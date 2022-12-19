import os, pickle, copy, copy
import pandas as pd
import numpy as np
import plotting.draw_biotypes as plots
from genomics.get_information_from_bed_intersection import GetInformationBEDIntersection
from collections import Counter
import operator

path_to_lib = '/'.join(os.path.abspath(__file__).split('/')[:-3])+'/lib/'

__author__ = "Maria Virginia Ruiz Cuevas"
__email__ = "maria.virginia.ruiz.cuevas@umontreal.ca"

class BiotypeAssignation:

	def __init__(self, path_to_output_folder, name_exp, mode, bam_files_list_rna, order_sample_bam_files_rna, dev, plots, super_logger, genome_version, mouse, threads):
		self.path_to_output_folder = path_to_output_folder
		self.name_exp = name_exp
		self.mode = mode
		self.plots = plots
		self.super_logger = super_logger
		self.mouse = mouse
		
		exists = os.path.exists(self.path_to_output_folder+'/res/BED_files/information_final_biotypes_peptides.dic') 
		exists_2 = os.path.exists(self.path_to_output_folder+'/res/BED_files/peptides_intersected_ere.dic') 
		
		if exists and exists_2:
				with open(self.path_to_output_folder+'/res/BED_files/peptides_intersected_ere.dic', 'rb') as handle:
					self.peptides_intersected_ere = pickle.load(handle)

				with open(self.path_to_output_folder+'/res/BED_files/information_final_biotypes_peptides.dic', 'rb') as handle:
					self.information_final_biotypes_peptides = pickle.load(handle)
		else:
			self.get_info_bed_files = GetInformationBEDIntersection(path_to_output_folder, self.mode, self.mouse, threads)
			self.get_info_bed_files.get_information_genomic_annotation(genome_version)
			self.get_info_bed_files.get_information_ERE_annotation()
			self.peptides_intersected_ere = self.get_info_bed_files.peptides_intersected_ere
			self.information_final_biotypes_peptides = self.get_info_bed_files.information_final_biotypes_peptides
			

		if not self.mouse:
			with open(path_to_lib+'ERE_info.dic', 'rb') as handle:
				self.ere_info = pickle.load(handle)
		else:
			with open(path_to_lib+'ERE_info_mouse.dic', 'rb') as handle:
				self.ere_info = pickle.load(handle)

		with open(path_to_lib+'coefficients.dic', 'rb') as handle:
			coefficients = pickle.load(handle)

		sorted_coefficients = dict( sorted(coefficients.items(), key=operator.itemgetter(1),reverse=True))
		self.biotypes_names = list(sorted_coefficients.keys())
		self.coefficients = list(sorted_coefficients.values())

		self.super_logger.info('========== Get information from Genomic and ERE annotation : Done! ============ ')

		self.order_sample_bam_files_rna = order_sample_bam_files_rna
		self.bam_files_list_rna = bam_files_list_rna
		self.dev = dev

			
	def get_biotypes(self, info_peptide_alignments, peptides_by_type):

		self.super_logger.info('========== Getting information to define biotyping... ============ ')
		data_gen_ere = []
		self.biotype_type = set()
		
		for type_peptide, peptides in peptides_by_type.items():

			for peptide in peptides:

				# EAAPDTVLR
				if True :#peptide == 'EAAPDTVLR'  : #== peptide :#'AEKLGFAGL' == peptide:
					try:
						info_alignments_peptide = info_peptide_alignments[peptide]
						alignments = info_alignments_peptide[0]
						total_count_rna = info_alignments_peptide[1]
						
						# ('chr13:99970439-99970465', 'GCGGCGGCGGCACCCCGGCCAGCTCTT', '+', [0, 0], [])
						# information_final_biotypes_peptides[key_peptide] = {transcript: [gene_type, transcript_type, transcript_level_biotype]}
						for key_peptide_alignment, info_alignment in alignments.items():

							peptide = key_peptide_alignment.split('_')[0]
							position = key_peptide_alignment.split('_')[1]
							MCS = key_peptide_alignment.split('_')[2]
							strand = info_alignment[0]
							rna_bam_files = info_alignment[1]
							
							to_add_gen_ere = []

							repName = ''
							repClass = ''
							repFamily = ''

							key_aux = peptide+'_'+position+'_'+MCS

							rep_names = []
							try:
								rep_names = list(self.peptides_intersected_ere[peptide][key_aux])
								repName = rep_names[0]
								if 'antisense_' in repName:
									repName_aux = repName.split('antisense_')[1]
									repClass = self.biotypes_names.index('Antisense_EREs')
									repFamily = self.ere_info[repName_aux][1]
									repFamily = 'antisense_'+repFamily
								else:
									repClass = self.biotypes_names.index(self.ere_info[repName][0])
									repFamily = self.ere_info[repName][1]
								
								self.biotype_type.add(repClass)

							except KeyError:
								pass

							try:
								transcripts_biotypes = self.information_final_biotypes_peptides[peptide][key_aux]

								for transcript, peptide_genomic_biotypes in transcripts_biotypes.items():
									
									gene_level_biotype = peptide_genomic_biotypes[0]
									transcript_level_biotype = peptide_genomic_biotypes[1]
									genomic_position_biotype = self.biotypes_names.index(self.mod_type(peptide_genomic_biotypes[2]))
										
									to_add_gen_ere = [type_peptide, peptide, position, MCS, strand, transcript, gene_level_biotype, transcript_level_biotype, genomic_position_biotype, repName, repClass, repFamily]
									
									to_add_gen_ere.extend(rna_bam_files)
									to_add_gen_ere.append(sum(rna_bam_files))
									data_gen_ere.append(to_add_gen_ere)
									
									self.biotype_type.add(genomic_position_biotype)

							except KeyError:
								transcript = 'No Annotation'
								gene_level_biotype = 'Intergenic'
								transcript_level_biotype = 'Intergenic'
								genomic_position_biotype = self.biotypes_names.index('Intergenic')

								if len(rep_names) > 0 :
									gene_level_biotype = ''
									transcript_level_biotype = ''
									genomic_position_biotype = ''
								else:
									self.biotype_type.add(genomic_position_biotype)
								
								to_add_gen_ere = [type_peptide, peptide, position, MCS, strand, transcript, gene_level_biotype, transcript_level_biotype, genomic_position_biotype, repName, repClass, repFamily]
									
								to_add_gen_ere.extend(rna_bam_files)
								to_add_gen_ere.append(sum(rna_bam_files))
								data_gen_ere.append(to_add_gen_ere)

					except KeyError:
						info_alignments_peptide = []

		columns_gen_ere = ['Peptide Type', 'Peptide','Alignment', 'MCS', 'Strand', 'Transcript', 'gene_level_biotype', 'transcript_level_biotype', 'genomic_position_biotype', 'ERE name', 'ERE class', 'ERE family']
		
		columns_gen_ere.extend(self.bam_files_list_rna)
		columns_gen_ere.extend(['Total reads count RNA'])
		self.data_gen_ere = pd.DataFrame(data_gen_ere, columns = columns_gen_ere)
		
		data_gen_ere = []
		self.super_logger.info('========== Getting information to define biotyping : Done! ============ ')


	def mod_type(self, type_):
		if '-' in type_:
			if 'Non_coding' in type_ :
				type_ = 'Non_coding Junctions'
			else:
				type_ = 'Junctions'
		return type_


	def get_global_annotation(self):

		self.bam_files_columns = self.bam_files_list_rna
		self.bam_files_columns.extend(['Total reads count RNA'])
		
		# Genomic and ERE Annotation Full - without MCS
		groupby_columns = ['Peptide Type', 'Peptide','Alignment', 'Strand', 'Transcript', 'gene_level_biotype', 'transcript_level_biotype', 'genomic_position_biotype', 'ERE name', 'ERE class', 'ERE family']
		
		self.df_total_by_position_gen_ere = self.data_gen_ere.groupby(groupby_columns)[self.bam_files_columns].sum().reset_index()
		self.df_total_by_position_gen_ere = pd.DataFrame(self.df_total_by_position_gen_ere, columns = groupby_columns+self.bam_files_columns)
		if self.dev:
			self.df_total_by_position_gen_ere.to_pickle(self.path_to_output_folder+'/res/biotype_classification/full_info_biotypes/df_total_by_position_gen_ere.pkl')
			
			
		def biotypes_translation(row):
			if row['genomic_position_biotype'] != '':
				row['genomic_position_biotype'] = self.biotypes_names[row['genomic_position_biotype']]
			if row['ERE class'] != '':
				if 'antisense_' in row['ERE name']:
					row['ERE class'] = 'antisense_'+self.biotypes_names[row['ERE class']]
				else:
					row['ERE class'] = self.biotypes_names[row['ERE class']]
			return row

		self.data_gen_ere = self.data_gen_ere.apply(lambda row : biotypes_translation(row), axis = 1)

		# Genomic and ERE Annotation Full - MCS
		path = self.path_to_output_folder+'/res/biotype_classification/full_info_biotypes/1_Genomic_and_ERE_Annotations_Full.csv'
		self.data_gen_ere.to_csv(path, index=False)
		
		del self.data_gen_ere

		df_total_by_position_gen_ere_to_print = self.df_total_by_position_gen_ere.apply(lambda row : biotypes_translation(row), axis = 1)
		path = self.path_to_output_folder+'/res/biotype_classification/full_info_biotypes/2_Genomic_and_ERE_Annotations_Summary_Full.csv'
		df_total_by_position_gen_ere_to_print.to_csv(path, index=False)
		
		del df_total_by_position_gen_ere_to_print

		self.super_logger.info('========== Global genomic and ERE annotation : Done! ============ ')
		
		self.get_genomic_and_transcription_based_biotype()
	

	def get_genomic_and_transcription_based_biotype(self):

		groupby_columns = ['Peptide Type',	'Peptide',	'Alignment', 'Strand']
		df_position = self.df_total_by_position_gen_ere.groupby(groupby_columns).count().reset_index()
		df_position = df_position.iloc[:, :len(groupby_columns)]
		self.df_position = pd.DataFrame(df_position, columns = groupby_columns)

		groupby_columns = ['Peptide Type',	'Peptide']
		df_position = self.df_total_by_position_gen_ere.groupby(groupby_columns).count().reset_index()
		df_position = df_position.iloc[:, :len(groupby_columns)]
		self.df_peptide = pd.DataFrame(df_position, columns = groupby_columns)
		
		groupby_columns = ['Peptide Type',	'Peptide',	'Alignment', 'Strand']
		df_position_genomic_position_ere_class = self.df_total_by_position_gen_ere.groupby(groupby_columns)
		
		groupby_columns = ['Peptide Type',	'Peptide',	'Alignment', 'Strand']+self.bam_files_columns
		counts = self.df_total_by_position_gen_ere.groupby(groupby_columns).count().reset_index()
		counts = counts.iloc[:, :len(groupby_columns)]
		self.counts = pd.DataFrame(counts, columns = groupby_columns)
		
		del self.df_total_by_position_gen_ere
		
		groupby_columns = ['Peptide Type',	'Peptide']
		total_count_by_alignment = counts.groupby(groupby_columns)[self.bam_files_columns].sum().reset_index()
		self.total_count_by_peptide = pd.DataFrame(total_count_by_alignment, columns = groupby_columns+self.bam_files_columns)
		
		df_position_biotypes_info_counts, df_position_biotypes_summary_genome = self.compute_genomic_and_ere_by_region(df_position_genomic_position_ere_class, counts)
		
		self.super_logger.info('========== Genomic and ERE annotation summary by position : Done! ============ ')
		print ('========== Genomic and ERE annotation summary by position : Done! ============ ')

		self.compute_genomic_and_ere_by_peptide(df_position_biotypes_summary_genome)

		self.super_logger.info('========== Genomic and ERE annotation by peptide : Done! ============ ')
		print ('========== Genomic and ERE annotation by peptide : Done! ============ ')

		extracted_col = self.compute_biotype_by_sample(df_position_biotypes_info_counts)

		self.super_logger.info('========== Transcription based biotyping by sample: Done! ============ ')
		print ('========== Transcription based biotyping by sample: Done! ============ ')

		self.compute_biotype_by_group_sample(df_position_biotypes_info_counts, extracted_col)

		self.super_logger.info('========== Transcription based biotyping by group sample: Done! ============ ')
		print ('========== Transcription based biotyping by group sample: Done! ============ ')


	def compute_genomic_and_ere_by_region(self, df_position_genomic_position_ere_class, counts):

		df_position_biotypes_info = copy.deepcopy(counts)
		df_position_biotypes_info.insert(4, 'Annotation Frequencies', '', allow_duplicates=False)
		df_position_biotypes_info.insert(5, 'Best Guess', '', allow_duplicates=False)

		df_position_biotypes_info_counts = pd.DataFrame(0.0, columns=range(len(self.biotypes_names)), index=range(len(self.df_position)))
		df_position_biotypes_summary_genome = copy.deepcopy(self.df_position)

		for index, biotype in enumerate(self.biotypes_names):
			df_position_biotypes_summary_genome[str(index)] = 0.0
			df_position_biotypes_summary_genome[str(index)] = pd.to_numeric(df_position_biotypes_summary_genome[str(index)], downcast="float")

		df_position_biotypes_summary_genome['Best Guess'] = ''

		columns_biotypes = df_position_biotypes_summary_genome.columns[2:-1]

	
		for index, row in self.df_position.iterrows():
			peptide_type = row[0]
			peptide = row[1]
			alignment = row[2]
			strand = row[3]

			group = df_position_genomic_position_ere_class.get_group((peptide_type, peptide, alignment, strand))

			genomic_biotypes = list(group['genomic_position_biotype'])
			ere_biotypes = list(filter(lambda a: a != '', group['ERE class']))
			total_biotypes = list(genomic_biotypes+ere_biotypes)


			if '' in total_biotypes:
				total_biotypes = list(filter(('').__ne__, total_biotypes))

			total_biotypes_types = len(total_biotypes)
			count_total_biotypes = Counter(total_biotypes)
			string_biotype = ''
			
			count_total_biotypes = dict(sorted(count_total_biotypes.items(), key = operator.itemgetter(1), reverse=True))
			
			best_guess = ''
			for bio, value in count_total_biotypes.items():
				
				coeff = self.coefficients[bio]
				df_position_biotypes_info_counts.at[index, bio] = coeff
				
				ratio = value/(total_biotypes_types*1.0)
				percentage = round(ratio*100,2)
				string_biotype += self.biotypes_names[bio]+': '+str(percentage)+'% - '
				
				df_position_biotypes_summary_genome.at[index, str(bio)] = ratio
				if self.biotypes_names[bio] == 'In_frame':
					best_guess = 'In_frame'

			string_biotype = string_biotype[:-2]
			if best_guess == '':
				best_guess = self.biotypes_names[list(count_total_biotypes.keys())[0]]#self.biotypes_names[sorted(total_biotypes)[0]]

			df_position_biotypes_info.at[index, 'Annotation Frequencies'] = string_biotype
			df_position_biotypes_info.at[index, 'Best Guess'] = best_guess

		groupby_columns = ['Peptide Type',	'Peptide']
		df_position_biotypes_summary_genome = df_position_biotypes_summary_genome.groupby(groupby_columns)[columns_biotypes].sum().reset_index()
		df_position_biotypes_summary_genome = pd.DataFrame(df_position_biotypes_summary_genome)

		path = self.path_to_output_folder+'/res/biotype_classification/full_info_biotypes/3_Genomic_and_ERE_Anno_by_Region_Full.csv'
		df_position_biotypes_info.to_csv(path, index=False)

		df_position_biotypes_info_counts = df_position_biotypes_info_counts.div(df_position_biotypes_info_counts.sum(axis=1), axis=0)
		
		return df_position_biotypes_info_counts, df_position_biotypes_summary_genome


	def compute_genomic_and_ere_by_peptide(self, df_position_biotypes_summary_genome):

		# Biotype based in the biotype genomic locations, no level of transcription!

		if self.plots:
			biotypes_by_peptide_type_genomic_ere_annot = {}
			biotypes_all_peptides_genomic_ere_annot = {}

		df_position_biotypes_info = copy.deepcopy(self.total_count_by_peptide)
		df_position_biotypes_info.insert(2, 'Annotation Frequencies', '', allow_duplicates=False)
		df_position_biotypes_info.insert(3, 'Best Guess', '', allow_duplicates=False)
		
		def get_biotype_without_transcription(row_input):
			
			peptide_type = row_input['Peptide Type']
			peptide = row_input['Peptide']
			
			row_summary = df_position_biotypes_summary_genome[(df_position_biotypes_summary_genome['Peptide Type'] == peptide_type) & (df_position_biotypes_summary_genome['Peptide'] == peptide)].values[0][2:-1]
			sum_alignments = sum(row_summary)
			
			final_values = (row_summary/sum_alignments)
			nonzeroind = np.nonzero(final_values)[0]
			b = sorted(enumerate(final_values), reverse=True, key=lambda i: i[1])[:len(nonzeroind)]
			
			string_biotype = ''
			best_guess = ''

			for index_biotype, ratio in b:
				percentage = round(ratio*100.0,2)
				biotype_name = self.biotypes_names[index_biotype]
				string_biotype += biotype_name+': '+str(percentage)+'% - '
				if biotype_name == 'In_frame':
					best_guess = biotype_name

				if self.plots:
					if ',' in peptide_type:
						types = peptide_type.split(',')
					else:
						types = peptide_type.split(';')

					for type_ in types:
						try:
							peptides_types = biotypes_by_peptide_type_genomic_ere_annot[type_]
							try:
								peptides_types[biotype_name] += ratio
							except KeyError:
								peptides_types[biotype_name] = ratio
						except KeyError:
							biotypes_by_peptide_type_genomic_ere_annot[type_] = {biotype_name:ratio}

					try:
						biotypes_all_peptides_genomic_ere_annot[biotype_name] += ratio
					except KeyError:
						biotypes_all_peptides_genomic_ere_annot[biotype_name] = ratio

			string_biotype = string_biotype[:-2]
			if best_guess == '' and len(b) > 0:
				best_guess = self.biotypes_names[b[0][0]]

			row_input['Annotation Frequencies'] = string_biotype
			row_input['Best Guess'] = best_guess

			df_position_biotypes_summary_genome.loc[(df_position_biotypes_summary_genome['Peptide Type'] == peptide_type) & (df_position_biotypes_summary_genome['Peptide'] == peptide), 'Best Guess'] = best_guess
			
			return row_input


		def get_biotype_without_transcription_summary(row_input):
			sum_alignments = sum(row_input[2:-1])
			row_input[2:-1] = (row_input[2:-1]/sum_alignments)
			return row_input

		df_position_biotypes_info = df_position_biotypes_info.apply(lambda row : get_biotype_without_transcription(row), axis = 1)
		path = self.path_to_output_folder+'/res/biotype_classification/summary_info_biotypes/1_General_Gen_and_ERE_Biotype_Consensus.csv'
		df_position_biotypes_info.to_csv(path, index=False)

		if self.dev:
			df_position_biotypes_summary_genome = df_position_biotypes_summary_genome.apply(lambda row : get_biotype_without_transcription_summary(row), axis = 1)
			df_position_biotypes_summary_genome.columns = ['Peptide Type', 'Peptide']+self.biotypes_names+['Best Guess']
		
			path = self.path_to_output_folder+'/res/biotype_classification/summary_info_biotypes/biotypes_by_peptide_genome_explained.csv'
			df_position_biotypes_summary_genome.to_csv(path, index=False)

		if self.plots:

			self.super_logger.info('========== Plots ============ ')
			
			plots.draw_biotypes(biotypes_by_peptide_type_genomic_ere_annot, self.path_to_output_folder+'plots/biotypes/genome_and_ERE_annotation/by_peptide_type/', False, False, self.name_exp)
			self.super_logger.info('========== biotypes_by_peptide_type_genomic_ere_annot : Genome & ERE annotations : Done! ============ ')
			biotypes_by_peptide_type_genomic_ere_annot = []

			plots.draw_biotypes(biotypes_all_peptides_genomic_ere_annot, self.path_to_output_folder+'plots/biotypes/genome_and_ERE_annotation/all_peptides/', True, False, self.name_exp)
			self.super_logger.info('========== biotypes_all_peptides_genomic_ere_annot : Genome & ERE annotations : Done! ============ ')
			biotypes_all_peptides_genomic_ere_annot = []

			self.super_logger.info('========== Plots : Done! ============ ')

	
	def compute_biotype_by_sample(self, df_position_biotypes_info_counts):

		if self.plots:
			biotypes_by_peptide_type_group_samples = {}
			biotypes_all_peptides_type_group_samples_all = {}

		df_biotype_by_peptide_by_sample = copy.deepcopy(self.df_peptide) 
		columns_biotypes = df_position_biotypes_info_counts.columns

		def get_biotype_with_transcription_by_sample(row_input):

			peptide = row_input['Peptide']
			
			try:
				total_count = int(self.total_count_by_peptide.loc[(self.total_count_by_peptide['Peptide'] == peptide)][bam_file])
			except TypeError:
				print (peptide, bam_file)
				exit()

			if total_count != 0:
				row_input[4:] = (row_input[4:]/total_count)
			else:
				row_input[4:] = 0

			return row_input

		def set_string_biotype_by_sample(row_input):

			peptide_type = row_input['Peptide Type']
			peptide = row_input['Peptide']

			row_summary = df_position_biotypes_info_2[(df_position_biotypes_info_2['Peptide Type'] == peptide_type) & (df_position_biotypes_info_2['Peptide'] == peptide)].values[0][2:-1]
			nonzeroind = np.nonzero(row_summary)[0]
			b = sorted(enumerate(row_summary), reverse=True, key=lambda i: i[1])[:len(nonzeroind)]
			
			string_biotype = ''
			best_guess = ''

			for index_biotype, ratio in b:
				percentage = round(ratio*100.0,2)
				biotype_name = self.biotypes_names[index_biotype]
				string_biotype += biotype_name+': '+str(percentage)+'% - '
				if biotype_name == 'In_frame':
					best_guess = biotype_name

				if self.plots:
					if ',' in peptide_type:
						types = peptide_type.split(',')
					else:
						types = peptide_type.split(';')
						

					if bam_file == 'Total reads count RNA':
						for type_ in types:
							try:
								all_samples = biotypes_by_peptide_type_group_samples['All_samples']
								try:
									peptides_types = all_samples[type_]
									try:
										peptides_types[biotype_name] += ratio
									except KeyError:
										peptides_types[biotype_name] = ratio
								except KeyError:
									all_samples[type_] = {biotype_name : ratio}
							except KeyError:
								dic_aux = {}
								dic_aux[type_] = {biotype_name : ratio}
								biotypes_by_peptide_type_group_samples['All_samples'] = dic_aux
						
						try:
							biotypes_all_peptides_type_group_samples_all[biotype_name] += ratio
						except KeyError:
							biotypes_all_peptides_type_group_samples_all[biotype_name] = ratio


			string_biotype = string_biotype[:-2]
			row_input[bam_file] = string_biotype
			
			if bam_file == 'Total reads count RNA':
				
				if best_guess == '' and len(b) > 0:
					best_guess = self.biotypes_names[b[0][0]]

				df_position_biotypes_info_2.loc[(df_position_biotypes_info_2['Peptide Type'] == peptide_type) & (df_position_biotypes_info_2['Peptide'] == peptide), 'Best Guess'] = best_guess
				
				if bam_file == 'Total reads count RNA' :
					row_input['Best Guess'] = best_guess
				
			return row_input

		for bam_file in ['Total reads count RNA']: #in self.bam_files_columns:
			df_biotype_by_peptide_by_sample[bam_file] = ''

			result = df_position_biotypes_info_counts.multiply(self.counts[bam_file], axis="index")
			result = pd.concat([self.df_position, result], axis=1)
			result = result.apply(lambda row : get_biotype_with_transcription_by_sample(row), axis = 1)
			
			groupby_columns = ['Peptide Type',	'Peptide']
			df_position_biotypes_info_2 = result.groupby(groupby_columns)[columns_biotypes].sum().reset_index()
			df_position_biotypes_info_2 = pd.DataFrame(df_position_biotypes_info_2)
			
			if bam_file == 'Total reads count RNA' :
				df_position_biotypes_info_2['Best Guess'] = ''
				df_biotype_by_peptide_by_sample['Best Guess'] = ''
			
			df_biotype_by_peptide_by_sample = df_biotype_by_peptide_by_sample.apply(lambda row : set_string_biotype_by_sample(row), axis = 1)
			
			if self.dev:
				if bam_file == 'Total reads count RNA':
					df_position_biotypes_info_2.columns = ['Peptide Type', 'Peptide']+self.biotypes_names+['Best Guess']
					path = self.path_to_output_folder+'/res/biotype_classification/summary_info_biotypes/biotypes_by_peptide_sample_explained_RNA.csv'
					df_position_biotypes_info_2.to_csv(path, index=False)

					peptide_alignment_sample = pd.concat([self.df_position, df_position_biotypes_info_counts, self.counts[bam_file]], axis=1)
					peptide_alignment_sample.drop(peptide_alignment_sample.index[peptide_alignment_sample['Total reads count RNA'] == 0], inplace = True)
					path = self.path_to_output_folder+'/res/biotype_classification/summary_info_biotypes/biotypes_by_peptide_alignment_and_sample_explained_RNA.csv'
					peptide_alignment_sample.to_csv(path, index=False)

		if self.plots:
			self.super_logger.info('========== Plots ============ ')

			plots.draw_biotypes(biotypes_by_peptide_type_group_samples, self.path_to_output_folder+'plots/biotypes/biotype_by_sample_group/by_peptide_type/', False, True, self.name_exp)
			self.super_logger.info('========== biotypes_by_peptide_type_group_samples : Transcription levels : Done! ============ ')
			biotypes_by_peptide_type_group_samples = []

			plots.draw_biotypes(biotypes_all_peptides_type_group_samples_all, self.path_to_output_folder+'plots/biotypes/biotype_by_sample_group/all_peptides/', True, False, self.name_exp)
			self.super_logger.info('========== biotypes_all_peptides_type_group_samples_all : Transcription levels : Done! ============ ')
			biotypes_all_peptides_type_group_samples_all = []

			self.super_logger.info('========== Plots : Done! ============ ')


		path = self.path_to_output_folder+'/res/biotype_classification/summary_info_biotypes/2_Sample_Gen_and_ERE_Biotype_Consensus.csv'
		df_biotype_by_peptide_by_sample.to_csv(path, index=False)
		extracted_cols = df_biotype_by_peptide_by_sample[['Total reads count RNA', 'Best Guess']]

		return extracted_cols
	

	def compute_biotype_by_group_sample(self, df_position_biotypes_info_counts, extracted_col):

		df_biotype_by_peptide_by_group_sample = copy.deepcopy(self.df_peptide) 
		columns_biotypes = df_position_biotypes_info_counts.columns

		total_count_by_peptide_by_group = copy.deepcopy(self.df_peptide)

		new_cols = ['Peptide Type', 'Peptide']+list(self.order_sample_bam_files_rna.keys())+['All']
		groups_columns = list(self.order_sample_bam_files_rna.keys())
		total_count_by_peptide_by_group = total_count_by_peptide_by_group.reindex(columns=new_cols, fill_value=0)
		df_biotype_by_peptide_by_group_sample = df_biotype_by_peptide_by_group_sample.reindex(columns=new_cols[:-1], fill_value='')
		
		def get_total_reads_by_group(row):
			peptide_type = row['Peptide Type']
			peptide = row['Peptide']

			col_names = list(row[2:].keys())

			for group in col_names:
				total_group = 0
				if group != 'All':
					samples = self.order_sample_bam_files_rna[group]
					for sample in samples:
						total_group += int(self.total_count_by_peptide.loc[(self.total_count_by_peptide['Peptide'] == peptide)][sample])
						
				row[group] = total_group

			row['All'] = int(self.total_count_by_peptide.loc[(self.total_count_by_peptide['Peptide'] == peptide)]['Total reads count RNA'])
			return row

		total_count_by_peptide_by_group = total_count_by_peptide_by_group.apply(lambda row : get_total_reads_by_group(row), axis = 1)
		
		def get_biotype_with_transcription_by_sample(row_input):

			peptide = row_input['Peptide']
			count_peptide_group = int(total_count_by_peptide_by_group.loc[(total_count_by_peptide_by_group['Peptide'] == peptide)][group])
			
			if count_peptide_group != 0:
				row_input[4:] = (row_input[4:]/count_peptide_group)	
			else:
				row_input[4:] = 0
			return row_input

		def set_string_biotype_by_group(row_input):

			peptide_type = row_input['Peptide Type']
			peptide = row_input['Peptide']

			row_summary = df_position_biotypes_info_2[(df_position_biotypes_info_2['Peptide Type'] == peptide_type) & (df_position_biotypes_info_2['Peptide'] == peptide)].values[0][2:-1]
			nonzeroind = np.nonzero(row_summary)[0]
			b = sorted(enumerate(row_summary), reverse=True, key=lambda i: i[1])[:len(nonzeroind)]
			
			string_biotype = ''
			
			for index_biotype, ratio in b:
				
				percentage = round(ratio*100.0,2)
				biotype_name = self.biotypes_names[index_biotype]
				string_biotype += biotype_name+': '+str(percentage)+'% - '
				

			string_biotype = string_biotype[:-2]
			row_input[group] = string_biotype

			return row_input

		for group in groups_columns:
			result = 0
			try:
				bam_files = self.order_sample_bam_files_rna[group]
			except KeyError:
				print (group, ' not in bamfiles ')

			for bam_file in bam_files:
				result += df_position_biotypes_info_counts.multiply(self.counts[bam_file], axis="index")
			
			result = pd.concat([self.df_position, result], axis=1)
			result = result.apply(lambda row : get_biotype_with_transcription_by_sample(row), axis = 1)
			
			groupby_columns = ['Peptide Type',	'Peptide']
			df_position_biotypes_info_2 = result.groupby(groupby_columns)[columns_biotypes].sum().reset_index()
			df_position_biotypes_info_2 = pd.DataFrame(df_position_biotypes_info_2)
			df_biotype_by_peptide_by_group_sample = df_biotype_by_peptide_by_group_sample.apply(lambda row : set_string_biotype_by_group(row), axis = 1)
			
		
		df_biotype_by_peptide_by_group_sample = df_biotype_by_peptide_by_group_sample.join(extracted_col)
		df_biotype_by_peptide_by_group_sample.rename(columns={"Total reads count RNA": "All"})
		
		path = self.path_to_output_folder+'/res/biotype_classification/summary_info_biotypes/3_Group_Samples_Gen_and_ERE_Biotype_Consensus.csv'
		df_biotype_by_peptide_by_group_sample.to_csv(path, index=False)

