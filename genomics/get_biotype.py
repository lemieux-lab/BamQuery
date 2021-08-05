import os, logging, time, pickle, multiprocessing, _thread, csv, math, copy, pysam, copy, xlsxwriter, gc
import pandas as pd
from pathos.multiprocessing import ProcessPool
import utils.useful_functions as uf
import plotting.plots as plots
from collections import Counter, OrderedDict
from genomics.get_information_from_bed_intersection import GetInformationBEDIntersection

path_to_lib = '/'.join(os.path.abspath(__file__).split('/')[:-3])+'/lib/'

__author__ = "Maria Virginia Ruiz Cuevas"
__email__ = "maria.virginia.ruiz.cuevas@umontreal.ca"

class BiotypeAssignation:

	def __init__(self, path_to_output_folder, name_exp, mode, bam_files_list_rna, bam_files_list_ribo, order_sample_bam_files_rna, order_sample_bam_files_ribo):
		self.path_to_output_folder = path_to_output_folder
		self.name_exp = name_exp
		self.mode = mode
		self.get_info_bed_files = GetInformationBEDIntersection(path_to_output_folder)
		self.get_info_bed_files.get_information_genomic_annotation()
		self.get_info_bed_files.get_information_ERE_annotation()

		with open(path_to_lib+'ERE_info.dic', 'rb') as handle:
			self.ere_info = pickle.load(handle)

		logging.info('========== Get information from Genomic and ERE annotation : Done! ============ ')

		self.order_sample_bam_files_rna = order_sample_bam_files_rna
		self.order_sample_bam_files_ribo = order_sample_bam_files_ribo
		self.bam_files_list_rna = bam_files_list_rna
		self.bam_files_list_ribo = bam_files_list_ribo
		
		
	def get_biotypes(self, info_peptide_alignments, peptides_by_type):

		logging.info('========== Getting information to define biotyping... ============ ')
		data_ere = []
		data_gen = []
		data_gen_ere = []
		self.biotype_type = set()

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
								self.biotype_type.add(repClass)

						except KeyError:
							pass

						try:
							transcripts_biotypes = self.get_info_bed_files.information_final_biotypes_peptides[peptide][key_aux]

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
								
								if '-' in genomic_position_biotype:
									if 'Non_coding' in genomic_position_biotype :
										self.biotype_type.add('Non_coding Junctions')
									else:
										self.biotype_type.add('Junctions')
								else:
									self.biotype_type.add(genomic_position_biotype)

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

		self.df1 = pd.DataFrame(data_ere, columns = columns_ere)
		self.df2 = pd.DataFrame(data_gen, columns = columns_gen)
		self.df3 = pd.DataFrame(data_gen_ere, columns = columns_gen_ere)
		data_gen_ere = []
		data_gen = []
		data_ere = []

		logging.info('========== Getting information to define biotyping : Done! ============ ')
		
	def get_global_annotation(self):

		# ERE Full
		groupby_columns = ['Peptide Type', 'Peptide','Alignment', 'Strand', 'ERE name', 'ERE class', 'ERE family']
		bam_files_columns = self.bam_files_list_rna
		bam_files_columns.extend(self.bam_files_list_ribo)
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
		
		#with open(self.path_to_output_folder+'alignments/df_total_by_position_gen_ere.dic', 'wb') as handle:
		#	pickle.dump(self.df_total_by_position_gen_ere, handle, protocol=pickle.HIGHEST_PROTOCOL)

		self.get_global_annotation_gen(bam_files_columns)
		#self.get_global_annotation_ere(bam_files_columns)
		#self.get_final_annotation_peptide_gen(bam_files_columns)
		#self.get_final_annotation_peptide_ere(bam_files_columns)

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

		#with open(self.path_to_output_folder+'alignments/df_global_gen.dic', 'wb') as handle:
		#	pickle.dump(self.df_global_gen, handle, protocol=pickle.HIGHEST_PROTOCOL)

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

		# with open(self.path_to_output_folder+'alignments/df_total_by_position_gen_ere.dic', 'rb') as handle:
		# 	self.df_total_by_position_gen_ere = pickle.load(handle)

		# with open(self.path_to_output_folder+'alignments/df_global_gen.dic', 'rb') as handle:
		# 	self.df_global_gen = pickle.load(handle)

		# print ('Charging !')
		# if bam_files_columns == '':
		# 	bam_files_columns = self.bam_files_list_rna
		# 	bam_files_columns.extend(self.bam_files_list_ribo)
		# 	bam_files_columns.extend(['Total reads count RNA', 'Total reads count Ribo'])

		logging.info('========== Init: Genomic and ERE annotation ============ ')

		df_consensus_annotation = []
		df_consensus_annotation_full = []
		df_consensus_annotation_final = []

		groupby_columns = ['Peptide Type', 'Peptide', 'Alignment', 'Strand']
		group_by_gen_ere = self.df_total_by_position_gen_ere.groupby(groupby_columns)

		groupby_columns = ['Peptide Type', 'Peptide']
		df_final_gen_sum = self.df_global_gen.groupby(groupby_columns)[bam_files_columns].sum()
		df_final_gen_sum = df_final_gen_sum.groupby(groupby_columns)

		peptides_visited = {}
		peptides_visited_alignment = {}
		peptides_visited_by_sample = {}
		peptides_visited_sample_group = {}

		biotypes_by_peptide_type = {}
		biotypes_all_peptides = {}
		
		biotypes_by_peptide_type_group_samples = {}
		biotypes_all_peptides_group_samples = {}
		biotypes_all_peptides_type_group_samples_all = {}
		
		counts_reads_rna_ribo = {}
		biotype_type = list(self.biotype_type)

		biotypes_by_peptide_genome_explained = []
		biotypes_by_peptide_sample_explained = []
		#biotype_type = ['Introns', 'Simple_repeat', 'In_frame', 'LINE', '3UTR', '5UTR', 'Non-coding Exons', 'Low_complexity', 'Intergenic', 'SINE', 'Junctions', 'Frameshift', 'LTR']

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

			if sum(count_bamfiles_alignment) > 0 :
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

			#if 'GSLGHLQNR' == peptide:
			alignment = key[2]
			count_bamfiles = df_final_gen_sum.get_group((peptide_type, peptide)).values.tolist()[0]
			count_bamfiles_alignment = group_by_gen_ere.get_group((key)).values.tolist()[0][11:]
			total_items_pep = 0
			
			type_gen_aux = list(filter(lambda a:a != '', items['genomic_position_biotype'].values.tolist()))
			type_ere_aux = list(filter(lambda a:a != '', items['ERE class'].values.tolist()))

			type_ere = []
			type_gen = []

			for type_ in type_ere_aux:
				types = type_.split(',')
				type_ere.extend(types)

			for type_ in type_gen_aux:
				types = type_.split(',')
				type_gen.extend(types)

			type_ere = list(set(type_ere))
			type_gen = list(set(type_gen))

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

			if sum(count_bamfiles_alignment) > 0 and len(peptides_visited_alignment[key]) == 1:
				peptides_visited_alignment[key].append(count_bamfiles_alignment)
				peptides_visited_alignment[key].append(count_bamfiles)
		

		aux_dic = {}
		biotype_info_aux_dic = {}

		peptides_absent_sample_group = {}
		
		#with open(self.path_to_output_folder+'alignments/peptides_visited.dic', 'wb') as handle:
		#	pickle.dump(peptides_visited, handle, protocol=pickle.HIGHEST_PROTOCOL)


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
			count_bamfiles_ = count_bamfiles[:-2]
			to_add_aux = []
			to_add_aux.extend(key)

			bios = [0]*len(biotype_type)
			
			if 0 in count_bamfiles_:
				indices = [i for i, x in enumerate(count_bamfiles_) if x == 0]
				peptides_absent_sample_group[peptide] = indices

			for type_ in sorted(dic, key=dic.get, reverse=True):
				count = dic[type_]
				count = count/(total_alignments_pep*1.0)
				percentage = round(count*100,2)
				count_string.append(type_+': '+str(percentage)+'%')
				
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


				if '-' in type_:
					if 'Non_coding' in type_ :
						type_ = 'Non_coding Junctions'
					else:
						type_ = 'Junctions'

				type_index = biotype_type.index(type_)
				bios[type_index] += percentage
				
				try:
					type_in_dic  = aux_dic['Genome']
					try:
						type_in_dic[type_]+= count
					except KeyError:
						type_in_dic[type_] = count
				except KeyError:
					aux_dic['Genome'] = {type_ : count}

				try:
					type_in_dic  = biotype_info_aux_dic['Genome']
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
					biotype_info_aux_dic['Genome'] = peptide_type_dic

			try:
				counts_reads_rna_ribo[peptide_type][peptide] = [count_rna, count_ribo]
			except KeyError:
				dic = {}
				dic[peptide] = [count_rna, count_ribo]
				counts_reads_rna_ribo[peptide_type] = dic
			
			to_add_aux.extend(bios)
			consensus = ' - '.join(count_string)
			to_add = [peptide_type, peptide, consensus]
			to_add.extend(count_bamfiles)
			df_consensus_annotation.append(to_add)
			biotypes_by_peptide_genome_explained.append(to_add_aux)

		with open(self.path_to_output_folder+'alignments/biotypes_by_peptide_genome_explained.list', 'wb') as handle:
			pickle.dump(biotypes_by_peptide_genome_explained, handle, protocol=pickle.HIGHEST_PROTOCOL)

		groupby_columns = ['Peptide Type', 'Peptide', 'Consenssus']
		self.df_consensus_annotation = pd.DataFrame(df_consensus_annotation, columns = groupby_columns+bam_files_columns)
		df_consensus_annotation = []

		split_bams_files, indexes_group, total_by_group, order_columns = self.get_info_group_samples(bam_files_columns)
		absent_group_peptide = {}
		
		for peptide, absents in peptides_absent_sample_group.items():
			for a in absents:
				key_group = indexes_group[a]
				try:
					pep = absent_group_peptide[key_group]
					try:
						pep[peptide] += 1
					except KeyError:
						pep[peptide] = 1
				except KeyError:
					dic_aux = {}
					dic_aux[peptide] = 1
					absent_group_peptide[key_group] = dic_aux


		for key, information_peptide in peptides_visited_alignment.items():
			
			to_add = []
			peptide_type = key[0]
			peptide = key[1]
			alignment = key[2]
			strand = key[3]
			total_alignments_pep = sum(peptides_visited_alignment[key][0].values())
			dic = peptides_visited_alignment[key][0]
			count_bamfiles = peptides_visited_alignment[key][1]
			count_bamfiles_total = peptides_visited_alignment[key][2]
			
			count_rna = count_bamfiles[-2]
			count_ribo = count_bamfiles[-1]

			count_string = []
			to_add = [peptide_type, peptide, alignment, strand]
			
			for type_ in sorted(dic, key=dic.get, reverse=True):
				count = dic[type_]
				count = count/(total_alignments_pep*1.0)
				percentage = str(round(count*100,2))
				count_string.append(type_+': '+percentage+'%')
				
				for index, bamfile in enumerate(count_bamfiles):
					value_in_bam_file = bamfile
					value_total_bam_file = count_bamfiles_total[index]
					key_group = indexes_group[index]
					total_by_key_group = total_by_group[key_group]

					try:
						value = count*(value_in_bam_file/value_total_bam_file)
					except ZeroDivisionError:
						value = 0

					try:
						absents = absent_group_peptide[key_group][peptide]
						try:
							value_group = value/((total_by_key_group - absents) * 1.0)
						except ZeroDivisionError:
							value_group = 0
					except KeyError:
						value_group = value/(total_by_key_group * 1.0)

					try:
						info_peptide = peptides_visited_by_sample[(peptide_type, peptide)]
						try:
							info_peptide[index][type_] += value
						except TypeError:
							dic_aux = {type_ : value}
							peptides_visited_by_sample[(peptide_type, peptide)][index] = dic_aux
						except KeyError:
							info_peptide[index][type_] = value
							
					except KeyError:
						dic_aux = [0]*len(count_bamfiles)
						dic_aux[index] = {type_ : value}
						peptides_visited_by_sample[(peptide_type, peptide)] = dic_aux
						
					try:
						info_peptide = peptides_visited_sample_group[(peptide_type, peptide)]
						try:
							info_peptide[key_group][type_] += value_group
						except TypeError:
							dic_aux = {type_ : value_group}
							peptides_visited_sample_group[(peptide_type, peptide)][key_group] = dic_aux
						except KeyError:
							info_peptide[key_group][type_] = value_group
					except KeyError:
						dic_aux = copy.deepcopy(split_bams_files)
						dic_aux[key_group] = {type_ : value_group}
						peptides_visited_sample_group[(peptide_type, peptide)] = dic_aux

					if key_group != 'Total reads count RNA' and key_group != 'Total reads count Ribo':
						
						try:
							types = biotypes_all_peptides_group_samples[key_group]
							try:
								types[type_] += value_group

							except KeyError:
								biotypes_all_peptides_group_samples[key_group][type_] = value_group
								
						except KeyError:
							biotypes_all_peptides_group_samples[key_group] = {type_ : value_group}

						try:
							peptide_types = biotypes_by_peptide_type_group_samples[key_group]
							try:
								types = peptide_types[peptide_type]
								try:
									types[type_] += value_group
								except KeyError:
									types[type_] = value_group
							except KeyError:
								peptide_types[peptide_type] = {type_ : value_group}
						except KeyError:
							dic_aux = {}
							dic_aux[peptide_type] = {type_ : value_group}
							biotypes_by_peptide_type_group_samples[key_group] = dic_aux

					if key_group == 'Total reads count RNA':
						try:
							biotypes_all_peptides_type_group_samples_all[type_] += value_group

						except KeyError:
							biotypes_all_peptides_type_group_samples_all[type_] = value_group

						if '-' in type_:
							if 'Non_coding' in type_ :
								type_ = 'Non_coding Junctions'
							else:
								type_ = 'Junctions'

						try:
							type_in_dic  = aux_dic['All']
							try:
								type_in_dic[type_]+= value
							except KeyError:
								type_in_dic[type_] = value
						except KeyError:
							aux_dic['All'] = {type_ : value}

						try:
							type_in_dic  = biotype_info_aux_dic['All']
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
							biotype_info_aux_dic['All'] = peptide_type_dic

			
			consensus = ' - '.join(count_string)
			to_add.append(consensus)
			to_add.extend(count_bamfiles)
			df_consensus_annotation_full.append(to_add)

		groupby_columns = ['Peptide Type', 'Peptide', 'Alignment', 'Strand', 'Consenssus']
		self.df_consensus_annotation_full = pd.DataFrame(df_consensus_annotation_full, columns = groupby_columns+bam_files_columns)

		self.write_xls_with_all_info_biotypes()
		del self.df3
		del self.df_total_by_position_gen_ere 
		del self.df_consensus_annotation_full 
		df_consensus_annotation_full = []
		

		for key, information_peptide in peptides_visited_by_sample.items():
			to_add_aux = []
			to_add_aux.extend(key)
			peptide_type = key[0]

			for index, bam_file_dic in enumerate(information_peptide):
				count_string = []
				to_add = []
				sample = bam_files_columns[index]
				to_add.append(sample)
				to_add.extend(key)
				bios = [0]*len(biotype_type)

				for type_ in sorted(bam_file_dic, key=bam_file_dic.get, reverse=True):
					count = bam_file_dic[type_]

					if '-' in type_:
						if 'Non_coding' in type_ :
							type_ = 'Non_coding Junctions'
						else:
							type_ = 'Junctions'

					type_index = biotype_type.index(type_)
					
					if count > 0:
						percentage = round(count*100,2)
						count_string.append(type_+': '+str(percentage)+'%')
						bios[type_index] += percentage

						if sample != 'Total reads count RNA':
							try:
								type_in_dic  = aux_dic[sample]
								try:
									type_in_dic[type_] += count
								except KeyError:
									type_in_dic[type_] = count
							except KeyError:
								aux_dic[sample] = {type_ : count}

							try:
								type_in_dic  = biotype_info_aux_dic[sample]
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
								biotype_info_aux_dic[sample] = peptide_type_dic


				to_add.extend(bios)
				consensus = ' - '.join(count_string)
				to_add_aux.append(consensus)
				biotypes_by_peptide_sample_explained.append(to_add)
			df_consensus_annotation_final.append(to_add_aux)

		with open(self.path_to_output_folder+'alignments/biotypes_by_peptide_sample_explained.list', 'wb') as handle:
			pickle.dump(biotypes_by_peptide_sample_explained, handle, protocol=pickle.HIGHEST_PROTOCOL)


		#with open(self.path_to_output_folder+'alignments/peptides_visited_sample_group.dic', 'wb') as handle:
		#	pickle.dump(peptides_visited_sample_group, handle, protocol=pickle.HIGHEST_PROTOCOL)

		with open(self.path_to_output_folder+'alignments/aux_dic.dic', 'wb') as handle:
			pickle.dump(aux_dic, handle, protocol=pickle.HIGHEST_PROTOCOL)

		with open(self.path_to_output_folder+'alignments/biotype_info_aux_dic.dic', 'wb') as handle:
			pickle.dump(biotype_info_aux_dic, handle, protocol=pickle.HIGHEST_PROTOCOL)
		
		self.get_biotype_by_group_sample(order_columns, peptides_visited_sample_group)

		groupby_columns = ['Peptide Type', 'Peptide']
		self.df_consensus_annotation_full_final = pd.DataFrame(df_consensus_annotation_final, columns = groupby_columns+bam_files_columns)
		df_consensus_annotation_final = []

		groupby_columns = ['Sample', 'Peptide Type', 'Peptide']
		self.biotypes_by_peptide_sample_explained = pd.DataFrame(biotypes_by_peptide_sample_explained, columns = groupby_columns+biotype_type)

		groupby_columns = ['Peptide Type', 'Peptide']
		self.biotypes_by_peptide_genome_explained = pd.DataFrame(biotypes_by_peptide_genome_explained, columns = groupby_columns+biotype_type)
		
		self.write_xls_info_biotypes_explained()

		self.write_xls_with_consensus_biotypes()
		logging.info('========== Writting out biotyping in xls files : Done! ============ ')

		del self.df_consensus_annotation 
		del self.df_consensus_annotation_full_final 
		del self.df_consensus_annotation_final_sample 
		gc.collect()

		logging.info('========== Plots ============ ')

		self.draw_biotypes(biotypes_by_peptide_type, self.path_to_output_folder+'plots/biotypes/genome_and_ERE_annotation/by_peptide_type/', False, False)
		logging.info('========== biotypes_by_peptide_type : Done! ============ ')
		self.draw_biotypes(biotypes_all_peptides, self.path_to_output_folder+'plots/biotypes/genome_and_ERE_annotation/all_peptides/', True, False)
		logging.info('========== biotypes_all_peptides : Done! ============ ')
		self.draw_biotypes(biotypes_all_peptides_type_group_samples_all, self.path_to_output_folder+'plots/biotypes/biotype_by_sample_group/all_peptides/', True, False)
		logging.info('========== biotypes_all_peptides_type_group_samples_all : Done! ============ ')
		self.draw_biotypes(biotypes_all_peptides_group_samples, self.path_to_output_folder+'plots/biotypes/biotype_by_sample_group/all_peptides/', True, True)
		logging.info('========== biotypes_all_peptides_group_samples : Done! ============ ')
		self.draw_biotypes(biotypes_by_peptide_type_group_samples, self.path_to_output_folder+'plots/biotypes/biotype_by_sample_group/by_peptide_type/', False, True)
		logging.info('========== biotypes_by_peptide_type_group_samples : Done! ============ ')
		
		logging.info('========== Plots : Done! ============ ')

		if self.mode == 'translation':
			self.draw_correlation(counts_reads_rna_ribo)

		logging.info('========== Fini: Genomic and ERE annotation : Done! ============ ')

	def get_info_group_samples(self, bam_files_columns):

		split_bams_files = {} 
		indexes_group = {}
		order_columns = []
		total_by_group = {}


		for group, samples in self.order_sample_bam_files_rna.items():
			order_columns.append(group)
			total_by_group[group] = len(samples)
			for sample in samples:
				index = bam_files_columns.index(sample)
				indexes_group[index] = group
				split_bams_files[group] = {}

		for group, samples in self.order_sample_bam_files_ribo.items():
			order_columns.append(group)
			total_by_group[group] = len(samples)
			for sample in samples:
				index = bam_files_columns.index(sample)
				indexes_group[index] = group
				split_bams_files[group] = {}

		indexes_group[len(bam_files_columns)-2] = 'Total reads count RNA'
		indexes_group[len(bam_files_columns)-1] = 'Total reads count Ribo'
		total_by_group['Total reads count RNA'] = 1
		total_by_group['Total reads count Ribo'] = 1

		split_bams_files['Total reads count RNA'] = {} 
		split_bams_files['Total reads count Ribo'] = {}
		order_columns.append('Total reads count RNA')
		order_columns.append('Total reads count Ribo')

		return split_bams_files, indexes_group, total_by_group, order_columns

	def get_biotype_by_group_sample(self, order_columns, peptides_visited_sample_group):

		df_consensus_annotation_final_sample = []
		#with open(self.path_to_output_folder+'alignments/peptides_visited_sample_group.dic', 'rb') as handle:
		#	peptides_visited_sample_group = pickle.load(handle)

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

	def draw_biotypes(self, biotypes_peptides, path_to_save, global_, samples):

		others_non_coding = ['retained_intron', 'bidirectional_promoter_lncRNA', 'transcribed_unitary_pseudogene', 'transcribed_unprocessed_pseudogene', 'sense_overlapping','processed_pseudogene', 'unprocessed_pseudogene']
		others_protein_coding = ['IG_V_gene', 'TEC', 'Exons']
		others_ere = ['DNA','RC', 'RNA','Satellite','Simple_repeat','Unknown', 'Low_complexity', 'rRNA','scRNA','snRNA','srpRNA','tRNA']
		

		organisation_labels = {'Protein-coding Regions':['5UTR', '3UTR', 'In_frame', 'Frameshift', 'protein_coding', 'CDS', 'Junctions', 'Other coding regions'], 
							'Non-coding Regions':['processed_transcript', 'nonsense_mediated_decay', 'antisense', 'Non_coding Exons', 'Non_coding Junctions', 'lincRNA', 'Other non_coding regions'], 
							'Intergenic Regions':['Intergenic'], 
							'Intronic Regions':['Introns'],
							'EREs':['LINE', 'LTR','Retroposon','SINE', 'Other EREs']}

		def plot_biotype(biotypes, name):
			title = name
			labels_in_type_peptide = {'Protein-coding genes':{}, 'Non-coding genes': {}, 'Protein-coding Regions':{}, 'Non-coding Regions': {}, 'Protein-coding transcripts':{}, 'Non-coding transcripts': {}, 'Intergenic Regions':{}, 'Intronic Regions':{}, 'EREs':{}}
			outer_labels = []
			outer_sizes = []
			intra_labels = []
			intra_sizes = []

			for biotype, total_biotype in biotypes.items():

				in_ = False
				for type_, types in organisation_labels.items():
					if biotype in types:
						labels_in_type_peptide[type_][biotype] = total_biotype 
						in_ = True
						break

				if not in_:
					
					if biotype in others_non_coding:
						type_ = 'Non-coding Regions'
						try:
							labels_in_type_peptide[type_]['Other non-coding regions'] += total_biotype 
						except KeyError:
							labels_in_type_peptide[type_]['Other non-coding regions'] = total_biotype

					elif biotype in others_protein_coding:
						type_ = 'Protein-coding Regions'
						try:
							labels_in_type_peptide[type_]['Other coding regions'] += total_biotype
						except KeyError:
							labels_in_type_peptide[type_]['Other coding regions'] = total_biotype

					elif biotype in others_ere:
						type_ = 'EREs'
						try:
							labels_in_type_peptide[type_]['Other EREs'] += total_biotype
						except KeyError:
							labels_in_type_peptide[type_]['Other EREs'] = total_biotype

					else:
						if '-' in biotype:
							if 'Non_coding' in biotype :
								try:
									labels_in_type_peptide['Non-coding Regions']['Non_coding Junctions'] += total_biotype 
								except KeyError:
									labels_in_type_peptide['Non-coding Regions']['Non_coding Junctions'] = total_biotype
							else:
								try:
									labels_in_type_peptide['Protein-coding Regions']['Junctions'] += total_biotype 
								except KeyError:
									labels_in_type_peptide['Protein-coding Regions']['Junctions'] = total_biotype
						else:
							print ('Problem with assignation biotype : ', biotype, name)
							
						
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
			gc.collect()

		if not samples:
			if not global_ :
				
				for type_peptide, biotypes in biotypes_peptides.items():
					plot_biotype(biotypes, type_peptide)
			else:
				plot_biotype(biotypes_peptides, 'All_peptides')
		else:
			if not global_  :
				for key_group, key_peptides_dic in biotypes_peptides.items():
					for type_peptide, biotypes in key_peptides_dic.items():
						plot_biotype(biotypes, type_peptide+'_'+key_group)
			else:
				for key_group, biotypes in biotypes_peptides.items():
					plot_biotype(biotypes, 'All_peptides_'+key_group)

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
			
		data = {'Type Peptide': peptide_type, 
				'Peptides': peptides,
        		'Read count Ribo': counts_ribo,
        		'Read count RNA': counts_rna} 
        		
		df_correlation = pd.DataFrame(data)
		plots.correlation(self.path_to_output_folder, self.name_exp, df_correlation)

	def write_xls_with_all_info_biotypes(self):

		writer = pd.ExcelWriter(self.path_to_output_folder+'/res/Annotation_Biotypes_full_info.xlsx', engine='xlsxwriter')

		writer.book.use_zip64()
		# self.df2.to_excel(writer, sheet_name='Genomic Annotation Biotypes')
		# self.df_total_by_position_gen.to_excel(writer, sheet_name='Genomic Annotation Biotypes_')
		# self.df_global_gen.to_excel(writer, sheet_name='Genomic Annotation By Region')
		# self.df_final_peptide_gen.to_excel(writer, sheet_name='Genomic Annotation By peptide')

		# worksheet1 = writer.sheets['Genomic Annotation Biotypes']
		# worksheet1.set_tab_color('green')
		# worksheet1 = writer.sheets['Genomic Annotation Biotypes_']
		# worksheet1.set_tab_color('green')
		# worksheet1 = writer.sheets['Genomic Annotation By Region']
		# worksheet1.set_tab_color('green')
		# worksheet1 = writer.sheets['Genomic Annotation By peptide']
		# worksheet1.set_tab_color('green')

		# self.df1.to_excel(writer, sheet_name='ERE Annotation Biotypes')
		# self.df_total_by_position_ere.to_excel(writer, sheet_name='ERE Annotation Biotypes_')
		# self.df_global_ere.to_excel(writer, sheet_name='ERE Annotation By Region')
		# self.df_final_peptide_ere.to_excel(writer, sheet_name='ERE Annotation By peptide')

		# worksheet1 = writer.sheets['ERE Annotation Biotypes']
		# worksheet1.set_tab_color('blue')
		# worksheet1 = writer.sheets['ERE Annotation Biotypes_']
		# worksheet1.set_tab_color('blue')
		# worksheet1 = writer.sheets['ERE Annotation By Region']
		# worksheet1.set_tab_color('blue')
		# worksheet1 = writer.sheets['ERE Annotation By peptide']
		# worksheet1.set_tab_color('blue')

		self.df3.to_excel(writer, sheet_name='Genomic & ERE Annotations')
		self.df_total_by_position_gen_ere.to_excel(writer, sheet_name='Genomic & ERE Annotations_')
		self.df_consensus_annotation_full.to_excel(writer, sheet_name='Genomic & ERE Anno. By Region')

		worksheet1 = writer.sheets['Genomic & ERE Annotations']
		worksheet1.set_tab_color('purple')
		worksheet1 = writer.sheets['Genomic & ERE Annotations_']
		worksheet1.set_tab_color('purple')
		worksheet1 = writer.sheets['Genomic & ERE Anno. By Region']
		worksheet1.set_tab_color('purple')

		#https://xlsxwriter.readthedocs.io/working_with_colors.html
		writer.save()

	def write_xls_info_biotypes_explained(self):
		
		writer = pd.ExcelWriter(self.path_to_output_folder+'/res/biotypes_by_peptide_sample_explained.xlsx', engine='xlsxwriter')

		writer.book.use_zip64()
		
		self.biotypes_by_peptide_sample_explained.to_excel(writer, sheet_name='Biotypes Sample Explained')
		self.biotypes_by_peptide_genome_explained.to_excel(writer, sheet_name='Biotypes Genome Explained')
		
		writer.save()

	def write_xls_with_consensus_biotypes(self):
		writer = pd.ExcelWriter(self.path_to_output_folder+'/res/Annotation_Biotypes_consensus.xlsx', engine='xlsxwriter')
		writer.book.use_zip64()

		self.df_consensus_annotation.to_excel(writer, sheet_name='General Gen & ERE Biotype')
		worksheet1 = writer.sheets['General Gen & ERE Biotype']
		worksheet1.set_tab_color('pink')

		self.df_consensus_annotation_full_final.to_excel(writer, sheet_name='Sample Gen & ERE Biotype')
		worksheet1 = writer.sheets['Sample Gen & ERE Biotype']
		worksheet1.set_tab_color('silver')
		
		self.df_consensus_annotation_final_sample.to_excel(writer, sheet_name='Group Samples Gen & ERE Biotype')
		worksheet1 = writer.sheets['Group Samples Gen & ERE Biotype']
		worksheet1.set_tab_color('navy')

		writer.save()


		
