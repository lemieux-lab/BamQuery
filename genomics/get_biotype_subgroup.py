import os, sys, time, pickle, _thread, csv, math, copy, pysam, copy, xlsxwriter, gc, argparse
import pandas as pd
import numpy as np
import billiard as mp
sys.path += sys.path + ['/u/ruizma/BAM_Query/Scripts/Python/BamQuery/plotting/']
import draw_biotypes as plots
from collections import Counter, OrderedDict
from collections import Counter
import operator

path_to_lib = '/'.join(os.path.abspath(__file__).split('/')[:-3])+'/lib/'

__author__ = "Maria Virginia Ruiz Cuevas"
__email__ = "maria.virginia.ruiz.cuevas@umontreal.ca"


class BiotypeAssignation:

	def __init__(self, path_to_output_folder, path_to_get_info_folder, name_exp, bam_files_list):
		self.path_to_output_folder = path_to_output_folder
		self.path_to_get_info_folder = path_to_get_info_folder
		self.name_exp = name_exp
		self.plots = plots
		
		output_paths = ['/full_info_biotypes/','/summary_info_biotypes/','/plots/biotypes/genome_and_ERE_annotation/by_peptide_type/','/plots/biotypes/genome_and_ERE_annotation/all_peptides/','/plots/biotypes/biotype_by_sample_group/by_peptide_type/','/plots/biotypes/biotype_by_sample_group/all_peptides/']

		for path in output_paths:
			try:
				os.makedirs(self.path_to_output_folder + path)
			except :
				pass

		with open(bam_files_list, 'rb') as fp:
			bam_files_info = pickle.load(fp)

		list_bam_files_order_rna = []
		order_sample_bam_files_rna = {}

		for name_sample, info_bam in sorted(bam_files_info.items(), key=lambda e: e[1][-2], reverse=False):
			list_bam_files_order_rna.append(name_sample)
			group = info_bam[-2]
			try:
				order_sample_bam_files_rna[group].append(name_sample)
			except KeyError:
				order_sample_bam_files_rna[group] = [name_sample]

		try:
			with open(self.path_to_get_info_folder+'/res/BED_files/peptides_intersected_ere.dic', 'rb') as handle:
				self.peptides_intersected_ere = pickle.load(handle)

			with open(self.path_to_get_info_folder+'/res/BED_files/information_final_biotypes_peptides.dic', 'rb') as handle:
				self.information_final_biotypes_peptides = pickle.load(handle)
		except ValueError:
			import pickle5
			with open(self.path_to_get_info_folder+'/res/BED_files/peptides_intersected_ere.dic', 'rb') as handle:
				self.peptides_intersected_ere = pickle5.load(handle)

			with open(self.path_to_get_info_folder+'/res/BED_files/information_final_biotypes_peptides.dic', 'rb') as handle:
				self.information_final_biotypes_peptides = pickle5.load(handle)

		with open(path_to_lib+'ERE_info.dic', 'rb') as handle:
			self.ere_info = pickle.load(handle)

		with open(path_to_lib+'coefficients.dic', 'rb') as handle:
			coefficients = pickle.load(handle)

		coefficients['Mutated'] = coefficients['Frameshift']
		sorted_coefficients = dict( sorted(coefficients.items(), key=operator.itemgetter(1),reverse=True))
		self.biotypes_names = list(sorted_coefficients.keys())
		self.coefficients = list(sorted_coefficients.values())

		with open(path_to_lib+'splices_information.dic', 'rb') as handle:
			self.splices_information = pickle.load(handle)

		self.order_sample_bam_files_rna = order_sample_bam_files_rna
		self.bam_files_list_rna = list_bam_files_order_rna
		self.bam_files_list_ribo = []
		self.order_sample_bam_files_ribo = {}
		self.dev = True


	def get_biotypes(self, info_peptide_alignments_path, peptides_by_type_path, peptides_to_filter):

		with open(peptides_to_filter, "rb") as fp:
			peptides_to_filter = pickle.load(fp)

		with open(info_peptide_alignments_path, 'rb') as fp:
			perfect_alignments = pickle.load(fp)

		with open(peptides_by_type_path, 'rb') as fp:
			peptides_by_type = pickle.load(fp)

		print ('Nb peptides to compute biotype classification : ', len(peptides_to_filter))

		def get_info_peptide_alignments():

			info_peptide_alignments = {}

			for peptide_alignment in perfect_alignments:
				peptide = peptide_alignment.split('_')[0]
				if peptide in peptides_to_filter:
					alignment = peptide_alignment.split('_')[1]
					MCS = peptide_alignment.split('_')[2]
					count_rna = perfect_alignments[peptide_alignment][-2]
					count_ribo = perfect_alignments[peptide_alignment][-1]
					strand = perfect_alignments[peptide_alignment][0]
					
					try:
						info_peptide = info_peptide_alignments[peptide]
						info_peptide[0][peptide_alignment] = [strand, count_rna, count_ribo]
						info_peptide[1] += sum(count_rna) 
						info_peptide[2] += sum(count_ribo)
					except KeyError:
						dic = {}
						dic[peptide_alignment] = [strand, count_rna, count_ribo]
						info_peptide_alignments[peptide] = [dic, sum(count_rna), sum(count_ribo)]

			return info_peptide_alignments

		info_peptide_alignments = get_info_peptide_alignments()

		data_gen_ere = []
		self.biotype_type = set()
		self.splices_annotated = set()

		for type_peptide, peptides in peptides_by_type.items():

			for peptide in peptides:

				if peptide in peptides_to_filter:
					try:
						info_alignments_peptide = info_peptide_alignments[peptide]

						alignments = info_alignments_peptide[0]
						total_count_rna = info_alignments_peptide[1]
						total_count_ribo = info_alignments_peptide[2]
						
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
									self.splices_annotated.add(position)

							to_add_gen_ere = []

							repNames = []
							repClass_s = []
							repFamilies = []

							repName = ''
							repClass = ''
							repFamily = ''

							key_aux = peptide+'_'+position

							rep_names = []
							try:
								rep_names = list(self.peptides_intersected_ere[peptide][key_aux])
								
								repName = rep_names[0]
								if 'antisense_' in repName:
									repName_aux = repName.split('antisense_')[1]
									repClass = self.biotypes_names.index('Antisense_'+self.biotypes_names[self.biotypes_names.index(self.ere_info[repName_aux][0])])
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
										
									#to_add_gen_ere = [type_peptide, peptide, position, MCS, strand, transcript, gene_level_biotype, transcript_level_biotype, genomic_position_biotype, ",".join(repNames), ",".join(repClass_s), ",".join(repFamilies)]
									to_add_gen_ere = [type_peptide, peptide, position, MCS, strand, transcript, gene_level_biotype, transcript_level_biotype, genomic_position_biotype, repName, repClass, repFamily]
									
									to_add_gen_ere.extend(rna_bam_files)
									to_add_gen_ere.extend(ribo_bam_files) 
									to_add_gen_ere.extend([sum(rna_bam_files), sum(ribo_bam_files)])
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
								
								#to_add_gen_ere = [type_peptide, peptide, position, MCS, strand, transcript, gene_level_biotype, transcript_level_biotype, genomic_position_biotype, ",".join(repNames), ",".join(repClass_s), ",".join(repFamilies)]
								to_add_gen_ere = [type_peptide, peptide, position, MCS, strand, transcript, gene_level_biotype, transcript_level_biotype, genomic_position_biotype, repName, repClass, repFamily]
									
								to_add_gen_ere.extend(rna_bam_files)
								to_add_gen_ere.extend(ribo_bam_files) 
								to_add_gen_ere.extend([sum(rna_bam_files), sum(ribo_bam_files)])
								data_gen_ere.append(to_add_gen_ere)

					except KeyError:
						info_alignments_peptide = []

		columns_gen_ere = ['Peptide Type', 'Peptide','Alignment', 'MCS', 'Strand', 'Transcript', 'gene_level_biotype', 'transcript_level_biotype', 'genomic_position_biotype', 'ERE name', 'ERE class', 'ERE family']
		columns_gen_ere.extend(self.bam_files_list_rna)
		columns_gen_ere.extend(self.bam_files_list_ribo)
		columns_gen_ere.extend(['Total reads count RNA', 'Total reads count Ribo'])

		self.data_gen_ere = pd.DataFrame(data_gen_ere, columns = columns_gen_ere)
		
		data_gen_ere = []
	
	
	def mod_type(self, type_):
		if '-' in type_:
			if 'Non_coding' in type_ :
				type_ = 'Non_coding Junctions'
			else:
				type_ = 'Junctions'
		return type_


	def get_global_annotation(self):

		self.bam_files_columns = self.bam_files_list_rna
		self.bam_files_columns.extend(self.bam_files_list_ribo)
		self.bam_files_columns.extend(['Total reads count RNA', 'Total reads count Ribo'])
		
		# Genomic and ERE Annotation Full - without MCS
		groupby_columns = ['Peptide Type', 'Peptide','Alignment', 'Strand', 'Transcript', 'gene_level_biotype', 'transcript_level_biotype', 'genomic_position_biotype', 'ERE name', 'ERE class', 'ERE family']
		self.df_total_by_position_gen_ere = self.data_gen_ere.groupby(groupby_columns)[self.bam_files_columns].sum().reset_index()
		self.df_total_by_position_gen_ere = pd.DataFrame(self.df_total_by_position_gen_ere, columns = groupby_columns+self.bam_files_columns)
		
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
		path = self.path_to_output_folder+'/full_info_biotypes/1_Genomic_and_ERE_Annotations_Full.csv'
		self.data_gen_ere.to_csv(path, index=False)
		
		del self.data_gen_ere

		df_total_by_position_gen_ere_to_print = self.df_total_by_position_gen_ere.apply(lambda row : biotypes_translation(row), axis = 1)
		path = self.path_to_output_folder+'/full_info_biotypes/2_Genomic_and_ERE_Annotations_Summary_Full.csv'
		df_total_by_position_gen_ere_to_print.to_csv(path, index=False)
		
		del df_total_by_position_gen_ere_to_print

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
		
		print ('========== Genomic and ERE annotation summary by position : Done! ============ ')

		self.compute_genomic_and_ere_by_peptide(df_position_biotypes_summary_genome)

		print ('========== Genomic and ERE annotation by peptide : Done! ============ ')

		extracted_col = self.compute_biotype_by_sample(df_position_biotypes_info_counts)

		print ('========== Transcription based biotyping by sample: Done! ============ ')

		self.compute_biotype_by_group_sample(df_position_biotypes_info_counts, extracted_col)

		print ('========== Transcription based biotyping by group sample: Done! ============ ')


	def compute_genomic_and_ere_by_region(self, df_position_genomic_position_ere_class, counts):

		df_position_biotypes_info = copy.deepcopy(counts)
		df_position_biotypes_info.insert(4, 'Annotation Frequencies', '', allow_duplicates=False)
		df_position_biotypes_info.insert(5, 'Best Guess', '', allow_duplicates=False)

		#df_position_biotypes_info_counts = copy.deepcopy(self.df_position)
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

		path = self.path_to_output_folder+'/full_info_biotypes/3_Genomic_and_ERE_Anno_by_Region_Full.csv'
		df_position_biotypes_info.to_csv(path, index=False)

		df_position_biotypes_info_counts = df_position_biotypes_info_counts.div(df_position_biotypes_info_counts.sum(axis=1), axis=0)
		
		return df_position_biotypes_info_counts, df_position_biotypes_summary_genome


	def compute_genomic_and_ere_by_peptide(self, df_position_biotypes_summary_genome):

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
			if best_guess == '':
				best_guess = self.biotypes_names[sorted(nonzeroind)[0]]

			row_input['Annotation Frequencies'] = string_biotype
			row_input['Best Guess'] = best_guess

			df_position_biotypes_summary_genome.loc[(df_position_biotypes_summary_genome['Peptide Type'] == peptide_type) & (df_position_biotypes_summary_genome['Peptide'] == peptide), 'Best Guess'] = best_guess
			
			return row_input


		def get_biotype_without_transcription_summary(row_input):
			sum_alignments = sum(row_input[2:-1])
			row_input[2:-1] = (row_input[2:-1]/sum_alignments)
			return row_input

		df_position_biotypes_info = df_position_biotypes_info.apply(lambda row : get_biotype_without_transcription(row), axis = 1)
		path = self.path_to_output_folder+'/summary_info_biotypes/1_General_Gen_and_ERE_Biotype_Consensus.csv'
		df_position_biotypes_info.to_csv(path, index=False)

		if self.dev:
			df_position_biotypes_summary_genome = df_position_biotypes_summary_genome.apply(lambda row : get_biotype_without_transcription_summary(row), axis = 1)
			df_position_biotypes_summary_genome.columns = ['Peptide Type', 'Peptide']+self.biotypes_names+['Best Guess']
		
			path = self.path_to_output_folder+'/summary_info_biotypes/biotypes_by_peptide_genome_explained.csv'
			df_position_biotypes_summary_genome.to_csv(path, index=False)

		if self.plots:

			plots.draw_biotypes(biotypes_by_peptide_type_genomic_ere_annot, self.path_to_output_folder+'/plots/biotypes/genome_and_ERE_annotation/by_peptide_type/', False, False, self.name_exp)
			biotypes_by_peptide_type_genomic_ere_annot = []
			
			plots.draw_biotypes(biotypes_all_peptides_genomic_ere_annot, self.path_to_output_folder+'/plots/biotypes/genome_and_ERE_annotation/all_peptides/', True, False, self.name_exp)
			biotypes_all_peptides_genomic_ere_annot = []

			
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
			
			if bam_file == 'Total reads count RNA' or bam_file == 'Total reads count Ribo':
				
				if best_guess == '' and len(b) > 0:
					best_guess = self.biotypes_names[b[0][0]]

				df_position_biotypes_info_2.loc[(df_position_biotypes_info_2['Peptide Type'] == peptide_type) & (df_position_biotypes_info_2['Peptide'] == peptide), 'Best Guess'] = best_guess
				
				if bam_file == 'Total reads count RNA' :
					row_input['Best Guess RNA'] = best_guess
				
				if bam_file == 'Total reads count Ribo':
					row_input['Best Guess Ribo'] = best_guess

			return row_input

		for bam_file in ['Total reads count RNA', 'Total reads count Ribo']: #in self.bam_files_columns:
			df_biotype_by_peptide_by_sample[bam_file] = ''

			result = df_position_biotypes_info_counts.multiply(self.counts[bam_file], axis="index")
			result = pd.concat([self.df_position, result], axis=1)
			result = result.apply(lambda row : get_biotype_with_transcription_by_sample(row), axis = 1)
			
			groupby_columns = ['Peptide Type',	'Peptide']
			df_position_biotypes_info_2 = result.groupby(groupby_columns)[columns_biotypes].sum().reset_index()
			df_position_biotypes_info_2 = pd.DataFrame(df_position_biotypes_info_2)
			
			if bam_file == 'Total reads count RNA' :
				df_position_biotypes_info_2['Best Guess'] = ''
				df_biotype_by_peptide_by_sample['Best Guess RNA'] = ''
			
			if bam_file == 'Total reads count Ribo':
				df_biotype_by_peptide_by_sample['Best Guess Ribo'] = ''

			df_biotype_by_peptide_by_sample = df_biotype_by_peptide_by_sample.apply(lambda row : set_string_biotype_by_sample(row), axis = 1)
			
			if self.dev:
				if bam_file == 'Total reads count RNA':
					df_position_biotypes_info_2.columns = ['Peptide Type', 'Peptide']+self.biotypes_names+['Best Guess']
					path = self.path_to_output_folder+'/summary_info_biotypes/biotypes_by_peptide_sample_explained_RNA.csv'
					df_position_biotypes_info_2.to_csv(path, index=False)

					peptide_alignment_sample = pd.concat([self.df_position, df_position_biotypes_info_counts, self.counts[bam_file]], axis=1)
					peptide_alignment_sample.drop(peptide_alignment_sample.index[peptide_alignment_sample['Total reads count RNA'] == 0], inplace = True)
					path = self.path_to_output_folder+'/summary_info_biotypes/biotypes_by_peptide_alignment_and_sample_explained_RNA.csv'
					peptide_alignment_sample.to_csv(path, index=False)

				if bam_file == 'Total reads count Ribo' :
					df_position_biotypes_info_2.columns = ['Peptide Type', 'Peptide']+self.biotypes_names+['Best Guess']
					path = self.path_to_output_folder+'/summary_info_biotypes/biotypes_by_peptide_sample_explained_Ribo.csv'
					df_position_biotypes_info_2.to_csv(path, index=False)


		if self.plots:
			
			plots.draw_biotypes(biotypes_by_peptide_type_group_samples, self.path_to_output_folder+'/plots/biotypes/biotype_by_sample_group/by_peptide_type/', False, True, self.name_exp)
			biotypes_by_peptide_type_group_samples = []

			plots.draw_biotypes(biotypes_all_peptides_type_group_samples_all, self.path_to_output_folder+'/plots/biotypes/biotype_by_sample_group/all_peptides/', True, False, self.name_exp)
			biotypes_all_peptides_type_group_samples_all = []


		path = self.path_to_output_folder+'/summary_info_biotypes/2_Sample_Gen_and_ERE_Biotype_Consensus.csv'
		df_biotype_by_peptide_by_sample.to_csv(path, index=False)
		extracted_cols = df_biotype_by_peptide_by_sample[['Total reads count RNA', 'Best Guess RNA', 'Total reads count Ribo', 'Best Guess Ribo']]

		return extracted_cols
	

	def compute_biotype_by_group_sample(self, df_position_biotypes_info_counts, extracted_col):

		df_biotype_by_peptide_by_group_sample = copy.deepcopy(self.df_peptide) 
		columns_biotypes = df_position_biotypes_info_counts.columns

		total_count_by_peptide_by_group = copy.deepcopy(self.df_peptide)

		new_cols = ['Peptide Type', 'Peptide']+list(self.order_sample_bam_files_rna.keys())+list(self.order_sample_bam_files_ribo.keys())+['All']
		groups_columns = list(self.order_sample_bam_files_rna.keys())+list(self.order_sample_bam_files_ribo.keys())
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
				bam_files = self.order_sample_bam_files_ribo[group]

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
		
		path = self.path_to_output_folder+'/summary_info_biotypes/3_Group_Samples_Gen_and_ERE_Biotype_Consensus.csv'
		df_biotype_by_peptide_by_group_sample.to_csv(path, index=False)


def main(argv):

	parser = argparse.ArgumentParser(description='======== Biotype Classification ========')
	
	parser.add_argument('path_to_output_folder', type=str, help='Path to the output folder where to save outputs information ')
	parser.add_argument('path_to_get_info_folder', type=str, help='Path to the output folder where to get the information : BEDs files ')
	parser.add_argument('name_exp', type=str, help='BamQuery search Id')
	parser.add_argument('bamfiles', type=str, help='bam_files_list')
	parser.add_argument('info_peptide_alignments_path', type=str, help='Get Information where the alignment information can be found ')
	parser.add_argument('peptides_by_type_path', type=str, help='Dic of all peptides in the original search : peptides_by_type_path')
	parser.add_argument('peptides_to_filter_path', type=str, help='Path of the list of peptides to be filtered for biotype computation analysis')
	
	args = parser.parse_args()
	
	path_to_output_folder = args.path_to_output_folder
	path_to_get_info_folder = args.path_to_get_info_folder
	name_exp = args.name_exp
	bam_files_list = args.bamfiles
	info_peptide_alignments_path = args.info_peptide_alignments_path
	peptides_by_type_path = args.peptides_by_type_path
	peptides_to_filter_path = args.peptides_to_filter_path


	t0 = time.time()

	biotype_classification = BiotypeAssignation(path_to_output_folder, path_to_get_info_folder, name_exp, bam_files_list)	
	biotype_classification.get_biotypes(info_peptide_alignments_path, peptides_by_type_path, peptides_to_filter_path)
	biotype_classification.get_global_annotation()

	t2 = time.time()
	total = t2-t0
	print ('Total time to run Biotype classification : ', str(total))	

if __name__ == "__main__":
	main(sys.argv[1:])


