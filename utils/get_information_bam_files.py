import os, threading, time, subprocess, concurrent.futures, getpass, pickle, sys, getopt, os, pysam, multiprocessing, csv
from os import listdir
from os.path import isfile, join
from pathos.multiprocessing import ProcessPool

__author__ = "Maria Virginia Ruiz"
__email__ = "maria.virginia.ruiz.cuevas@umontreal.ca"

path_to_lib = '/'.join(os.path.abspath(__file__).split('/')[:-3])+'/lib/'

files_with_not_permission = []

NUM_WORKERS =  int(multiprocessing.cpu_count()/2)

class GetInformationBamFiles:

	def __init__(self, path_to_input_folder, path_to_output_folder, mode, strandedness, light, bam_files_logger):
		
		path_to_lib = '/'.join(os.path.abspath(__file__).split('/')[:-3])+'/lib/'
		path_to_all_counts_file = path_to_lib+"allcounts.dic"
		exist = os.path.exists(path_to_all_counts_file)
		
		self.bam_files_list = {}
		self.bam_ribo_files_list = {}
		self.bam_files_logger = bam_files_logger

		if light:
			self.path_to_output_temps_folder = path_to_output_folder+'res_light/temps_files/'
			self.path_to_output_aux_folder = path_to_output_folder+'res_light/AUX_files/'
		else:
			self.path_to_output_temps_folder = path_to_output_folder+'res/temps_files/'
			self.path_to_output_aux_folder = path_to_output_folder+'res/AUX_files/'
		
		
		try:
			bam_files = path_to_input_folder+'BAM_directories.tsv'
			self.bam_files_list = self.get_info_bamfiles(bam_files, strandedness, path_to_output_folder)
		
		except FileNotFoundError:
			self.bam_files_logger.info('The bam directories : %s doesn\'t exist ', path_to_input_folder+'BAM_directories.tsv')

		if mode == 'translation':
			try:
				bam_files = path_to_input_folder+'BAM_Ribo_directories.tsv'
				self.bam_ribo_files_list = self.get_info_bamfiles(bam_files, strandedness, path_to_output_folder)

			except FileNotFoundError:
				self.bam_files_logger.info('If running translation mode you must include a list of Ribo Bam Files. The bam directories : %s doesn\'t exist ', path_to_input_folder+'BAM_Ribo_directories.tsv')


	def get_info_bamfiles(self, bam_files, strandedness, path_to_output_folder):

		data = [['Sample', 'Tissue', 'Tissue_type', 'Short_list']]
		
		bam_files_list = {}
		exists = os.path.exists(self.path_to_output_temps_folder+"bam_files_info_query.dic")
		initial_list_paths = []

		mod = False
		bam_files_to_get_primary_read_count = []

		if not exists:

			path_to_all_counts_file = path_to_lib+"Bam_files_info.dic"
			exist = os.path.exists(path_to_all_counts_file)
			
			if exist :
				try:
					with open(path_to_all_counts_file, 'rb') as fp:
						dictionary_total_reads_bam_files = pickle.load(fp)
				except ValueError:
					import pickle5
					with open(path_to_all_counts_file, 'rb') as fp:
						dictionary_total_reads_bam_files = pickle5.load(fp)
			else:
				dictionary_total_reads_bam_files = {}
				print ('Bam_files_info dictionary is not found in the lib path. If this is the first time you are running \
					BAMQuery or you are querying new samples, this will take a bit of time. Go for coffee! If it is not \
					the case you should be worried, the Bam_files_info dictionary have been lost, now is generating a new Bam_files_info dictionary \
					with the samples in this query.')
				self.bam_files_logger.info('Bam_files_info dictionary is not found in the lib path. If this is the first time you are running \
					BAMQuery or you are querying new samples, this will take a bit of time. Go for coffee! If it is not \
					the case you should be worried, the Bam_files_info dictionary have been lost, now is generating a new Bam_files_info dictionary \
					with the samples in this query.')

			
			with open(bam_files) as f:

				for index, line in enumerate(f):
					line = line.strip().split('\t')
					name = line[0]

					try:
						path = line[1]
						initial_list_paths.append(path)
					except IndexError:
						raise Exception("Sorry, your BAM_directories.tsv file does not follow the correct format. Remember that the columns must be tab separated.")


					if '.bam' not in path or '.cram' not in path:
						bam_files_found = self.search_bam_files(path)
					else:
						bam_files_found = [path]


					for bam_file_path in bam_files_found:

						name_bam_file = "_".join(bam_file_path.split('/')[:-1][-2:])
						sequencing = ''
						library = ''
						# 0: path, 1: number, 2: Tissue, 3: Tissue_type, 4: Shortlist, 5: sequencing, 6: library, 7: user
						try:
							info_bam_file = dictionary_total_reads_bam_files[name_bam_file]
							path = info_bam_file[0]
							count = info_bam_file[1]
							tissue = info_bam_file[2]
							tissue_type = info_bam_file[3]
							shortlist = info_bam_file[4]
							sequencing = info_bam_file[5]
							library = info_bam_file[6]
							user = info_bam_file[7]
							if count == 0:
								bam_files_to_get_primary_read_count.append(bam_file_path)
							if tissue == '' or tissue_type == '' or shortlist == '':
								to_add = [''] * 4
								to_add[0] = name_bam_file
								data.append(to_add)
						except KeyError:
							path = bam_file_path
							count = 0
							tissue = ''
							tissue_type = ''
							shortlist = ''
							sequencing = ''
							library = ''
							user = getpass.getuser()
							info_bam_file = [bam_file_path, count, tissue, tissue_type, shortlist, sequencing, library, user]
							dictionary_total_reads_bam_files[name_bam_file] = info_bam_file

							to_add = [''] * 4
							to_add[0] = name_bam_file
							data.append(to_add)

							bam_files_to_get_primary_read_count.append(bam_file_path)
							mod = True
							
						if strandedness:
							if sequencing == '' and library == '':
								library, sequencing = self.get_type_library(bam_file_path)
								dictionary_total_reads_bam_files[name_bam_file][5] = sequencing
								dictionary_total_reads_bam_files[name_bam_file][6] = library
								mod = True	
						else:
							library = 'unstranded' 
							sequencing = 'unstranded' 

						if library == '' and sequencing == '':
							self.bam_files_logger.info('Sample bam file for %s %s is not properly indexed, skipping...', name, name_bam_file) 
						else:
							bam_files_list[name_bam_file] = [bam_file_path, sequencing, library, name]

						if path != bam_file_path:
							self.bam_files_logger.info('Information Change: BAM file is already in the dictionary, however the path for the BAM file is not the same. \
							Path in BAM_directories.tsv: %s, Path already assigned %s. Total reads : %s', bam_file_path, path, str(count))
							dictionary_total_reads_bam_files[name_bam_file][0] = path
							mod = True
						else :
							self.bam_files_logger.info('%s total reads : %s ', bam_file_path, str(count))


			with open(self.path_to_output_temps_folder+"bam_files_info_query.dic", 'wb') as handle:
				pickle.dump(bam_files_list, handle, protocol=pickle.HIGHEST_PROTOCOL)

			if mod :
				with open(path_to_all_counts_file, 'wb') as handle:
					pickle.dump(dictionary_total_reads_bam_files, handle, protocol=pickle.HIGHEST_PROTOCOL)

			for bam_file in files_with_not_permission:
				if bam_file in initial_list_paths:
					self.bam_files_logger.info('Skipping BAM file : %s. This Bam file has access denied or the path does not exist..', bam_file) 

		else:
			
			try:
				with open(self.path_to_output_temps_folder+"bam_files_info_query.dic", 'rb') as fp:
					bam_files_list = pickle.load(fp)
			except ValueError:
				import pickle5
				with open(self.path_to_output_temps_folder+"bam_files_info_query.dic", 'rb') as fp:
					bam_files_list = pickle5.load(fp)


		self.bam_files_logger.info('Total Bam Files to Query : %d.', len(bam_files_list))

		if len(bam_files_to_get_primary_read_count) > 0 :
				
			path_to_save_bam_files_to_search = self.path_to_output_aux_folder+"bam_files_to_get_primary_read_count.dic"
			
			with open(path_to_save_bam_files_to_search, 'wb') as handle:
				pickle.dump(bam_files_to_get_primary_read_count, handle, protocol=pickle.HIGHEST_PROTOCOL)

			get_read_counts_path = '/'.join(os.path.abspath(__file__).split('/')[:-1])+'/primary_read_count.py'

			command = 'python '+get_read_counts_path+' -i '+path_to_save_bam_files_to_search+' -o '+ path_to_output_folder
			subprocess.Popen(command, shell=True, stdout=subprocess.DEVNULL, close_fds=True)
			self.bam_files_logger.info('Total Bam Files to get primary read counts : %d ', len(data) - 1)

		if len(data) > 1:
			self.bam_files_logger.info('Please enter the tissue information for the new BamFiles into the %s file. ', self.path_to_output_aux_folder+'bam_files_tissues.csv')
			
			with open(self.path_to_output_aux_folder+"bam_files_tissues.csv", 'w') as csvFile:
				writer = csv.writer(csvFile)
				writer.writerows(data)

			raise Exception("\nBefore to continue you must provide the tissue type for the bam files annotated in the file : "+ self.path_to_output_aux_folder+"bam_files_tissues.csv. Please enter for each sample : tissue, tissue_type, shortlist." )

		return bam_files_list
	

	def search_bam_files(self, path):

		bam_files_found = []

		for p, d, f in os.walk(path, onerror=self.walk_error_handler, followlinks=True):

			for file in f:
				if (file.endswith('.cram') and (file+'.crai' in f or file[:-5]+'.bai' in f)) or (file.endswith('.bam') and (file+'.bai' in f or file[:-4]+'.bai' in f)):
					bam_file = join(p, file)
					bam_files_found.append(bam_file)

		return bam_files_found

	def walk_error_handler(self, exception_instance):
		files_with_not_permission.append(exception_instance.filename)


	def get_type_library(self, path):
		
		try:
			count_1 = int(pysam.view("-f1",'-c', path, 'chr12:6,537,097-6,537,227'))

			if count_1 == 0 :
				sequencing = 'single-end'
				count_1 = int(pysam.view("-f0X10",'-c', path, 'chr12:6,537,097-6,537,227'))
				count_2 = int(pysam.view("-F0X10",'-c', path, 'chr12:6,537,097-6,537,227'))
			else:
				sequencing = 'pair-end'
				count_1 = int(pysam.view("-f0X50",'-c', path, 'chr12:6,537,097-6,537,227')) # Conversion -f80 to hexa
				count_2 = int(pysam.view("-f0X60",'-c', path, 'chr12:6,537,097-6,537,227')) # Conversion -f96 to hexa

			ratio = 0
			try:
				ratio = (count_1+count_2)/(abs(count_1-count_2)*1.0)
			except :
				pass
			
			type_library = ''
			
			if ratio > 2 :
				type_library = 'unstranded'
			else:
				if count_1 > count_2:
					type_library = 'reverse'
				elif count_2 > count_1:
					type_library = 'forward'
				elif count_1 == count_2 and count_1 == 0:
					self.bam_files_logger.info('Guessing library for this Bam file %s fail. Adding unstranded library ! ' , path)
					type_library = 'unstranded'

			return type_library, sequencing

		except pysam.utils.SamtoolsError: 

			return '',''



