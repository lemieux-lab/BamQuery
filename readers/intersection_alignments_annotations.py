import os, time, subprocess


__author__ = "Maria Virginia Ruiz Cuevas"

class IntersectAnnotations:

	def __init__(self, perfect_alignments, path_to_output_folder, mode, name_exp, super_logger, genome_version, mouse):
		path_to_lib = '/'.join(os.path.abspath(__file__).split('/')[:-3])+'/lib/'
		self.perfect_alignments = perfect_alignments
		self.mouse = mouse

		if mode == 'translation':
			self.path_to_output_folder = path_to_output_folder+'res_translation/BED_files/'
		else:
			self.path_to_output_folder = path_to_output_folder+'res/BED_files/'
		
		self.bed_file = self.path_to_output_folder + 'to_intersect_to_annotations.bed'
		self.name_exp = name_exp
		self.super_logger = super_logger
		
		if self.mouse:
			self.annotation_EREs = path_to_lib + 'EREs_souris.bed'
			if genome_version == 'M24':
				self.annotation_transcripts = path_to_lib+'genome_versions/genome_mouse_m24/gencode.vM24.primary_assembly.annotation.gtf'
			if genome_version == 'M30':
				self.annotation_transcripts = path_to_lib+'genome_versions/genome_mouse_m30/gencode.vM30.primary_assembly.annotation.gtf'
			
		else:
			self.annotation_EREs = path_to_lib + 'hg38_ucsc_repeatmasker.gtf'

			if genome_version == 'v26_88': 
				self.annotation_transcripts = path_to_lib+'genome_versions/genome_v26_88/gencode.v26.primary_assembly.annotation.gtf'
			elif genome_version == 'v33_99':
				self.annotation_transcripts = path_to_lib+'genome_versions/genome_v33_99/gencode.v33.primary_assembly.annotation.gtf'
			elif genome_version == 'v38_104':
				self.annotation_transcripts = path_to_lib+'genome_versions/genome_v38_104/gencode.v38.primary_assembly.annotation.gtf'

		

	def generate_BED_files(self):

		exists = os.path.exists(self.path_to_output_folder+ 'to_intersect_to_annotations.bed')
		
		if not exists:
			self.super_logger.info('Generating BED file to intersect with annotations. ')

			t_0 = time.time()

			to_write = ''
			places_written = set()

			for key, info_alignment in self.perfect_alignments.items():
				split_key = key.split('_')
				peptide = split_key[0]
				chr = split_key[1].split(':')[0]
				places = split_key[1].split(':')[1]
				try:
					MCS = split_key[2]
				except IndexError:
					pass

				strand = info_alignment[0]
				key_unique = peptide+'_'+split_key[1]

				if key_unique not in places_written:

					if '|' in places:
						split_places = places.split('|')
						for index, place in enumerate(split_places): # chr1:28214940-28214963
							to_write += self.get_string_for_alignment(place, chr, key_unique, strand)
					else:
						to_write += self.get_string_for_alignment(places, chr, key_unique, strand)

					places_written.add(key_unique)

			
			file_to_save = open(self.bed_file, 'w')
			file_to_save.write(to_write)
			file_to_save.close()

			t_2 = time.time()
			total = t_2-t_0
			self.super_logger.info('Total time run function generate_BED_files to end : %f min', (total/60.0))

		else:
			self.super_logger.info('to_intersect_to_annotations.bed file already collected in the output folder : %s --> Skipping this step!', self.path_to_output_folder + 'to_intersect_to_annotations.bed')
			
	
	def perform_intersection_with_annotation(self):

		exists = os.path.exists(self.path_to_output_folder+ 'intersection_with_annotated_transcripts.bed')
		exists_2 = os.path.exists(self.path_to_output_folder+ 'intersection_with_annotated_EREs.bed')
		
		if not exists and not exists_2:
			t_0 = time.time()

			self.super_logger.info('Using bedtools to intersect alignments to annotations.')

			command = 'bedtools intersect -a '+self.bed_file+' -b '+self.annotation_transcripts+' -wao | grep -w transcript > '+self.path_to_output_folder + '/intersection_with_annotated_transcripts.bed; bedtools intersect -a '+self.bed_file+' -b '+self.annotation_EREs+' -wao > '+self.path_to_output_folder + '/intersection_with_annotated_EREs.bed'
			p_1 = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
			out, err = p_1.communicate()

			t_2 = time.time()
			total = t_2-t_0
			self.super_logger.info('Total time run function perform_intersection_with_annotation to end : %f min', (total/60.0))
		else:
			self.super_logger.info('intersection_with_annotations.bed file already collected in the output folder : %s --> Skipping this step!', self.path_to_output_folder + 'intersection_with_annotations.bed')
	
	def perform_intersection_with_annotation_mouse(self):

		exists = os.path.exists(self.path_to_output_folder+ 'intersection_with_annotated_transcripts.bed')
		
		if not exists :
			t_0 = time.time()

			self.super_logger.info('Using bedtools to intersect alignments to annotations.')

			command = 'bedtools intersect -a '+self.bed_file+' -b '+self.annotation_transcripts+' -wao | grep -w transcript > '+self.path_to_output_folder + '/intersection_with_annotated_transcripts.bed; bedtools intersect -a '+self.bed_file+' -b '+self.annotation_EREs+' -wao > '+self.path_to_output_folder + '/intersection_with_annotated_EREs.bed' 
			p_1 = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
			out, err = p_1.communicate()

			t_2 = time.time()
			total = t_2-t_0
			self.super_logger.info('Total time run function perform_intersection_with_annotation to end : %f min', (total/60.0))
		else:
			self.super_logger.info('intersection_with_annotations.bed file already collected in the output folder : %s --> Skipping this step!', self.path_to_output_folder + 'intersection_with_annotations.bed')
	

	def perform_intersection_with_assemblies(self, assemblies_to_intersect):

		t_0 = time.time()

		self.super_logger.info('Using bedtools to intersect alignments to assemblies.')

		command = ''
		for assembly, bed_to_intersect in assemblies_to_intersect:
			command += 'bedtools intersect -a '+self.bed_file+' -b '+assembly+' -wao | grep -w StringTie > '+bed_to_intersect + ';'
		p_1 = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
		out, err = p_1.communicate()

		t_2 = time.time()
		total = t_2-t_0
		self.super_logger.info('Total time run function perform_intersection_with_annotation to end : %f min', (total/60.0))


	def get_string_for_alignment(self, place, chr, key, strand):
		start = place.split('-')[0]
		end = place.split('-')[1]
		string = chr+'\t'+start+'\t'+end+'\t'+key+'\t'+strand+'\n'
		
		return string

