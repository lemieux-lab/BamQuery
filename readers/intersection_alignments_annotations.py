import os, logging, time, pickle, multiprocessing, _thread, csv, math, subprocess
import pandas as pd

path_to_lib = '/'.join(os.path.abspath(__file__).split('/')[:-3])+'/lib/'

annotation_transcripts = path_to_lib + 'gencode.v26.primary_assembly.annotation.gtf'
annotation_EREs = path_to_lib + 'hg38_ucsc_repeatmasker.gtf'

__author__ = "Maria Virginia Ruiz Cuevas"

class IntersectAnnotations:

	def __init__(self, perfect_alignments, path_to_output_folder, name_exp):
		path_to_lib = '/'.join(os.path.abspath(__file__).split('/')[:-3])+'/lib/'
		self.perfect_alignments = perfect_alignments
		self.path_to_output_folder = path_to_output_folder+'res/BED_files/'
		self.name_exp = name_exp

	def generate_BED_files(self):

		exists = os.path.exists(self.path_to_output_folder+ 'to_intersect_to_annotations.bed')
		
		if not exists:
			logging.info('Generating BED file to intersect with annotations. ')

			t_0 = time.time()

			to_write = ''
			cont = 0
			for key, info_alignment in self.perfect_alignments.items():
				split_key = key.split('_')
				peptide = split_key[0]
				chr = split_key[1].split(':')[0]
				places = split_key[1].split(':')[1]
				strand = info_alignment[0]
				#TTRPALQEL_chrX:155355661-155355687 ['+', 'ACCACCAGACCTGCCCTACAAGAGCTC', 'TTRPALQEL', []]

				if '|' in places:
					split_places = places.split('|')
					for index, place in enumerate(split_places): # chr1:28214940-28214963
						id = key+'_'+str(cont)+'_'+str(index)
						to_write += self.get_string_for_alignment(place, chr, key, strand)
				else:
					id = key+'_'+str(cont)+'_0'
					to_write += self.get_string_for_alignment(places, chr, key, strand)

				cont += 1

			self.bed_file = self.path_to_output_folder + 'to_intersect_to_annotations.bed'
			file_to_save = open(self.bed_file, 'w')
			file_to_save.write(to_write)
			file_to_save.close()

			t_2 = time.time()
			total = t_2-t_0
			logging.info('Total time run function generate_BED_files to end : %f min', (total/60.0))

		else:
			logging.info('to_intersect_to_annotations.bed file already collected in the output folder : %s --> Skipping this step!', self.path_to_output_folder + 'to_intersect_to_annotations.bed')
			self.bed_file = self.path_to_output_folder + 'to_intersect_to_annotations.bed'

	def perform_intersection_with_annotation(self):

		exists = os.path.exists(self.path_to_output_folder+ 'intersection_with_annotated_transcripts.bed')
		exists_2 = os.path.exists(self.path_to_output_folder+ 'intersection_with_annotated_EREs.bed')
		
		if not exists and not exists_2:
			t_0 = time.time()

			logging.info('Using bedtools to intersect alignments to annotations.')

			command = 'module add bedtools; bedtools intersect -a '+self.bed_file+' -b '+annotation_transcripts+' -wao | grep -w transcript > '+self.path_to_output_folder + '/intersection_with_annotated_transcripts.bed; bedtools intersect -a '+self.bed_file+' -b '+annotation_EREs+' -wao > '+self.path_to_output_folder + '/intersection_with_annotated_EREs.bed'
			p_1 = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
			out, err = p_1.communicate()

			t_2 = time.time()
			total = t_2-t_0
			logging.info('Total time run function perform_intersection_with_annotation to end : %f min', (total/60.0))
		else:
			logging.info('intersection_with_annotations.bed file already collected in the output folder : %s --> Skipping this step!', self.path_to_output_folder + 'intersection_with_annotations.bed')
	

	def get_string_for_alignment(self, place, chr, key, strand):
		start = place.split('-')[0]
		end = place.split('-')[1]
		string = chr+'\t'+start+'\t'+end+'\t'+key+'\t'+strand+'\n'
		return string

