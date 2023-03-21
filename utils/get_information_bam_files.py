import os, subprocess, getpass, pickle, csv
from os.path import join
import pandas as pd
import inspect, sys, pysam
import time

__author__ = "Maria Virginia Ruiz"
__email__ = "maria.virginia.ruiz.cuevas@umontreal.ca"

path_to_lib = '/'.join(os.path.abspath(__file__).split('/')[:-3])+'/lib/'

files_with_not_permission = []


class NoTraceBackWithLineNumber(Exception):
    def __init__(self, msg):
        try:
            ln = sys.exc_info()[-1].tb_lineno
        except AttributeError:
            ln = inspect.currentframe().f_back.f_lineno
        self.args = "{0.__name__} (line {1}): {2}".format(type(self), ln, msg),
        sys.exit(self)

class NeedMoreInfo(NoTraceBackWithLineNumber):
    pass


class GetInformationBamFiles:

	def __init__(self, path_to_input_folder, path_to_output_folder, mode, strandedness, light, bam_files_logger, sc, genome_version, mouse, threads):
		
		self.bam_files_list = {}
		self.bam_ribo_files_list = {}
		self.bam_files_logger = bam_files_logger
		self.sc = sc
		self.mouse = mouse
		self.threads = threads
		self.path_to_input_folder = path_to_input_folder
		self.path_to_lock_file = path_to_lib+"lock_dic"
			

		if genome_version == 'v26_88': 
			self.genome_version_gtf= path_to_lib+'genome_versions/genome_v26_88/gencode.v26.primary_assembly.annotation.gtf'
		elif genome_version == 'v33_99':
			self.genome_version_gtf= path_to_lib+'genome_versions/genome_v33_99/gencode.v33.primary_assembly.annotation.gtf'
		else:
			self.genome_version_gtf = path_to_lib+'genome_versions/genome_v38_104/gencode.v38.primary_assembly.annotation.gtf'
		

		if light:
			self.path_to_output_temps_folder = path_to_output_folder+'res_light/temps_files/'
			self.path_to_output_aux_folder = path_to_output_folder+'res_light/AUX_files/'
		else:
			self.path_to_output_temps_folder = path_to_output_folder+'res/temps_files/'
			self.path_to_output_aux_folder = path_to_output_folder+'res/AUX_files/'
		
		if mode == 'translation':
			self.path_to_output_temps_folder = path_to_output_folder+'res_translation/temps_files/'
			self.path_to_output_bed_files_folder = path_to_output_folder+'res_translation/BED_files/'
			self.path_to_output_aux_folder = path_to_output_folder+'res_translation/AUX_files/'

			bam_files = path_to_input_folder+'BAM_Ribo_directories.tsv'
			exists = os.path.exists(bam_files)
			if exists:
				self.bam_files_list = self.get_info_bamfiles(bam_files, strandedness, path_to_output_folder)
			else:
				self.bam_files_logger.info('If running translation mode you must include a list of Ribo Bam Files. The BAM_directories.tsv file %s doesn\'t exist, please provide this file and relaunch the query.', bam_files)
		else:
			bam_files = path_to_input_folder+'BAM_directories.tsv'
			exists = os.path.exists(bam_files)
			if exists:
				self.bam_files_list = self.get_info_bamfiles(bam_files, strandedness, path_to_output_folder)
			else:
				self.bam_files_logger.info('The BAM_directories.tsv file %s doesn\'t exist, please provide this file and relaunch the query.', bam_files)
				message = '\nThe BAM_directories.tsv file '+ bam_files+' doesn\'t exist, please provide this file and relaunch the query.'
				raise NeedMoreInfo(message)


	def get_info_ribo_bamfiles(self, bam_files):
		
		ribo_bam_files_info_query = self.path_to_output_temps_folder+"ribo_bam_files_info_query.dic"
		exists = os.path.exists(ribo_bam_files_info_query)
		
		if not exists:
			bam_files_list = {}
			command = 'module load stringtie/1.3.6; ulimit -s 8192;'

			with open(bam_files) as f:

				for index, line in enumerate(f):
					line = line.strip().split('\t')
					name = line[0]

					try:
						path = line[1]
					except IndexError:
						raise Exception("Your BAM_directories.tsv file does not follow the correct format. Remember that the columns must be tab separated.")

					if '.bam' in path or '.cram' in path:
						bam_files_found = [path]
					else:
						bam_files_found = self.search_bam_files(path)
					
					assemblage = 0
					for bam_file_path in bam_files_found:

						name_bam_file = "_".join(bam_file_path.split('/')[:-1][-2:])
						output_assembled_gtf = self.path_to_output_temps_folder+name_bam_file+'.gtf'
						output_assembled_bed = self.path_to_output_bed_files_folder+name_bam_file+'.bed'
						bam_files_list[name_bam_file] = [bam_file_path, output_assembled_gtf, output_assembled_bed]

						exists_gtf = os.path.exists(output_assembled_gtf)
						
						if not exists_gtf:
							assemblage += 1 
							command += 'stringtie '+bam_file_path+' -G '+self.genome_version_gtf+' -o '+output_assembled_gtf+' -c 1 -p 32 --fr -m 30 -g 28 ;'
			
			with open(ribo_bam_files_info_query, 'wb') as handle:
				pickle.dump(bam_files_list, handle, protocol=pickle.HIGHEST_PROTOCOL)

			if assemblage > 0:
				subprocess.Popen(command, shell=True, stdout=subprocess.DEVNULL, close_fds=True).wait()
				self.bam_files_logger.info('Total Bam Files to get transcriptome assembly : %d ', assemblage)

		else:
			with open(ribo_bam_files_info_query, 'rb') as fp:
				bam_files_list = pickle.load(fp)

			assemblage = 0
			command = 'module load stringtie/1.3.6; ulimit -s 8192;'
			for name_bam_file, bam_file_info in bam_files_list.items():
				output_assembled_gtf = bam_file_info[1]
				bam_file_path = bam_file_info[0]
				exists = os.path.exists(output_assembled_gtf)
				if not exists:
					command += 'stringtie '+bam_file_path+' -G '+self.genome_version_gtf+' -o '+output_assembled_gtf+' -c 1 -p 32 --fr -m 30 -g 28;'
					assemblage += 1 

			if assemblage > 0:
				subprocess.Popen(command, shell=True, stdout=subprocess.DEVNULL, close_fds=True).wait()
				self.bam_files_logger.info('Total Bam Files to get transcriptome assembly : %d ', assemblage)

		return bam_files_list


	def remove_lock_to_bam_files_info_dic(self):
		with open(self.path_to_lock_file, 'r') as file:
			first_line = file.readline().strip()
		if first_line == self.path_to_input_folder:
			first_line = ''
			with open(self.path_to_lock_file, 'w') as file:
				file.writelines(first_line) 
			file.close()
		else:
			print ('Error in the path_to_lock_file: contents ', first_line, ' and the ', self.path_to_input_folder, ' not the same.')
		

	def grant_access_to_bam_files_info_dic(self, timeout=300):
		start_time = time.time()
		with open(self.path_to_lock_file, 'r') as file:
			first_line = file.readline().strip()
		while (time.time() - start_time < timeout) and first_line != self.path_to_input_folder:
			try:
				with open(self.path_to_lock_file, 'r') as file:
					first_line = file.readline().strip()
				if first_line == '':
					with open(self.path_to_lock_file, 'w') as file:
						file.write(self.path_to_input_folder)
					file.close()
					return True
				time.sleep(1)
			except IOError:
				pass
		raise NeedMoreInfo("\nTimeout to grant access to Bam_files_info.dic. Another BamQuery process may be using the dictionary." )


	def get_info_bamfiles(self, bam_files, strandedness, path_to_output_folder):

		data = [['Sample', 'Tissue', 'Tissue_type', 'Short_list']]
		
		bam_files_list = {}
		exists = os.path.exists(self.path_to_output_temps_folder+"bam_files_info_query.dic")
		initial_list_paths = []
		exists_lock_file = os.path.exists(self.path_to_lock_file)

		# First time creating the file lock_dic file
		if not exists_lock_file:
			file_to_open = open(self.path_to_lock_file, 'w')
			file_to_open.write('')
			file_to_open.close()
			self.bam_files_logger.info('Creating a lock file for managed access to Bam_files_info')

		mod = False
		bam_files_to_get_primary_read_count = []

		if not exists:

			path_to_all_counts_file = path_to_lib+"Bam_files_info.dic"
			exist = os.path.exists(path_to_all_counts_file)
			
			if self.grant_access_to_bam_files_info_dic():
				self.bam_files_logger.info('Lock access to Bam_files_info.dic')
				pass
			else:
				raise NeedMoreInfo("\nCould not acquire the lock to consult the Bam_files_info.dic. Check if there is another run of BamQuery that you wrote to the lock file and it is still running." )
			
			if exist :
				with open(path_to_all_counts_file, 'rb') as fp:
					dictionary_total_reads_bam_files = pickle.load(fp)
			else:
				dictionary_total_reads_bam_files = {}
				print ('The Bam_files_info dictionary is not in the library path. If this is the first time you are running BamQuery or you are querying new samples, this will take a little time while the primary read count for each BAM file is retrieved. If this is not the case, the Bam_files_info dictionary has been lost and BamQuery is now generating a new Bam_files_info dictionary with the samples from this query.')
				self.bam_files_logger.info('The Bam_files_info dictionary is not in the library path. If this is the first time you are running BamQuery or you are querying new samples, this will take a little time while the primary read count for each BAM file is retrieved. If this is not the case, the Bam_files_info dictionary has been lost and BamQuery is now generating a new Bam_files_info dictionary with the samples from this query.')

			
			with open(bam_files) as f:

				for index, line in enumerate(f):
					line = line.strip().split('\t')
					name = line[0]

					try:
						path = line[1]
						initial_list_paths.append(path)
					except IndexError:
						self.remove_lock_to_bam_files_info_dic()
						self.bam_files_logger.info('Unlock Bam_files_info')

						raise Exception("Your BAM_directories.tsv file does not follow the correct format. Remember that the columns must be tab separated.")

					if '.bam' in path or '.cram' in path:
						bam_files_found = [path]
					else:
						bam_files_found = self.search_bam_files(path)

					self.bam_files_logger.info('Total number of accessible bam files in the path %s  : %s ', path, str(len(bam_files_found)))

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
							library = '' 
							sequencing = ''
							user = getpass.getuser()
							info_bam_file = [bam_file_path, count, tissue, tissue_type, shortlist, sequencing, library, user]
							dictionary_total_reads_bam_files[name_bam_file] = info_bam_file

							to_add = [''] * 4
							to_add[0] = name_bam_file
							data.append(to_add)

							mod = True
						
						if strandedness:
							if (sequencing == '' and library == '') or (sequencing == 'unstranded' and library == 'unstranded'):
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
							bam_files_list[name_bam_file] = [bam_file_path, sequencing, library, name, count, tissue, tissue_type]

						if path != bam_file_path:
							self.bam_files_logger.info('Information Change: BAM file is already in the dictionary, however the path for the BAM file is not the same. Path in BAM_directories.tsv: %s, Path already assigned %s. Total reads : %s', path, bam_file_path, str(count))
							dictionary_total_reads_bam_files[name_bam_file][0] = bam_file_path
							mod = True
						else :
							self.bam_files_logger.info('%s total reads : %s ', bam_file_path, str(count))
					self.bam_files_logger.info('-----------')

			with open(self.path_to_output_temps_folder+"bam_files_info_query.dic", 'wb') as handle:
				pickle.dump(bam_files_list, handle, protocol=pickle.HIGHEST_PROTOCOL)

			if mod :
				with open(path_to_all_counts_file, 'wb') as handle:
					pickle.dump(dictionary_total_reads_bam_files, handle, protocol=pickle.HIGHEST_PROTOCOL)

			self.remove_lock_to_bam_files_info_dic()
			self.bam_files_logger.info('Unlock Bam_files_info')

			in_ = False	
			for bam_file in files_with_not_permission:
				if bam_file in initial_list_paths:
					self.bam_files_logger.info('Skipping BAM file : %s. This Bam file has access denied or the path does not exist..', bam_file) 
					in_ = True
			if in_:
				raise NeedMoreInfo("\nBefore proceeding, please verify that you have the access granted to query all BAM files or that the path to these files exists. \nSee the log file Information_BAM_directories.log for more information about BAM files with limited access." )
		else:
			with open(self.path_to_output_temps_folder+"bam_files_info_query.dic", 'rb') as fp:
				bam_files_list = pickle.load(fp)

			# 0: path, 1: number, 2: Tissue, 3: Tissue_type, 4: Shortlist, 5: sequencing, 6: library, 7: user
			if os.path.exists(self.path_to_output_aux_folder+"bam_files_tissues.csv"):
				
				if self.grant_access_to_bam_files_info_dic():
					self.bam_files_logger.info('Lock access to Bam_files_info.dic')
					pass
				else:
					raise NeedMoreInfo("\nCould not acquire the lock to consult the Bam_files_info.dic. Check if there is another run of BamQuery that you wrote to the lock file and it is still running." )
			
				self.bam_files_logger.info('Adding tissue information to Bam Files !')
				path_to_all_counts_file = path_to_lib+"Bam_files_info.dic"
			
				with open(path_to_all_counts_file, 'rb') as fp:
					dictionary_total_reads_bam_files = pickle.load(fp)

				df = pd.read_csv(self.path_to_output_aux_folder+"bam_files_tissues.csv", header = 0)

				for index, row in df.iterrows():
					try:
						sample = str(row['Sample'])
						tissue_name = str(row['Tissue'])
						tissue_type = str(row['Tissue_type'])
						short_list = str(row['Short_list']).lower()
						if sample != 'nan' and tissue_name != 'nan' and tissue_type != 'nan' and (short_list != 'yes' or short_list != 'no'):
							dictionary_total_reads_bam_files[sample][2] = tissue_name
							dictionary_total_reads_bam_files[sample][3] = tissue_type
							dictionary_total_reads_bam_files[sample][4] = short_list
							bam_files_list[sample].append(tissue_name)
							bam_files_list[sample].append(tissue_type)
						else:
							self.remove_lock_to_bam_files_info_dic()
							self.bam_files_logger.info('Unlock Bam_files_info')
							raise NeedMoreInfo("\nBefore to continue you must provide all the tissue information for the bam files annotated in the file : "+ self.path_to_output_aux_folder+"bam_files_tissues.csv. Please enter for each sample : tissue, tissue_type, shortlist (yes/no)." )
					except Exception :
						self.remove_lock_to_bam_files_info_dic()
						self.bam_files_logger.info('Unlock Bam_files_info')
						raise NeedMoreInfo("\nBefore to continue you must provide all the tissue information for the bam files annotated in the file : "+ self.path_to_output_aux_folder+"bam_files_tissues.csv. Please enter for each sample : tissue, tissue_type, shortlist (yes/no)." )
				
				with open(path_to_all_counts_file, 'wb') as handle:
					pickle.dump(dictionary_total_reads_bam_files, handle, protocol=pickle.HIGHEST_PROTOCOL)

				self.remove_lock_to_bam_files_info_dic()
				self.bam_files_logger.info('Unlock Bam_files_info')

				with open(self.path_to_output_temps_folder+"bam_files_info_query.dic", 'wb') as handle:
					pickle.dump(bam_files_list, handle, protocol=pickle.HIGHEST_PROTOCOL)

			for sample, info_sample in bam_files_list.items():
				bam_file_path = info_sample[0]
				count = info_sample[4]
				if count == 0:
					bam_files_to_get_primary_read_count.append(bam_file_path)

			self.bam_files_logger.info('Total Bam Files to Query : %d.', len(bam_files_list))
			
			if len(bam_files_to_get_primary_read_count) > 0 and not self.sc:
				
				path_to_save_bam_files_to_search = self.path_to_output_aux_folder+"bam_files_to_get_primary_read_count.dic"
				
				with open(path_to_save_bam_files_to_search, 'wb') as handle:
					pickle.dump(bam_files_to_get_primary_read_count, handle, protocol=pickle.HIGHEST_PROTOCOL)

				get_read_counts_path = '/'.join(os.path.abspath(__file__).split('/')[:-1])+'/primary_read_count.py'
				command = 'python '+get_read_counts_path+' -i '+path_to_save_bam_files_to_search+' -o '+ path_to_output_folder+' -t '+str(self.threads)
				subprocess.Popen(command, shell=True, stdout=subprocess.DEVNULL, close_fds=True)
				self.bam_files_logger.info('Total Bam Files to get primary read counts : %d ', len(bam_files_to_get_primary_read_count))

		if len(data) > 1 and not self.sc:
			self.bam_files_logger.info('Please enter the tissue information for the new BamFiles into the %s file. ', self.path_to_output_aux_folder+'bam_files_tissues.csv')
			
			with open(self.path_to_output_aux_folder+"bam_files_tissues.csv", 'w') as csvFile:
				writer = csv.writer(csvFile)
				writer.writerows(data)

			raise NeedMoreInfo("\nBefore to continue you must provide all tissue information for the bam files annotated in the file : "+ self.path_to_output_aux_folder+"bam_files_tissues.csv. Please enter for each sample : tissue, tissue_type, shortlist (yes/no)." )

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

		# MOUSE GAPDH Chromosome 6: 125,138,678-125,143,430 
		
		if self.mouse:
			position_gapdh_grch38 = 'chr6:125,163,139-125,163,436'
			position_gapdh_grch39 = 'chr6:125,138,678-125,143,430'
			try:
				count_gapdh_grch38 = int(pysam.view("-f1",'-c', path, position_gapdh_grch38))
				count_gapdh_grch39 = int(pysam.view("-f1",'-c', path, position_gapdh_grch39))
			except pysam.utils.SamtoolsError: 
				return '',''
			
			if count_gapdh_grch38 > count_gapdh_grch39:
				print ('GRCH 38')
				sequencing = 'pair-end'
				count_2 = int(pysam.view("-f0X50",'-c', path, position_gapdh_grch38)) # Conversion -f80 to hexa
				count_1 = int(pysam.view("-f0X60",'-c', path, position_gapdh_grch38)) # Conversion -f96 to hexa
				print ('-f80 ',count_2, '-f96 ',count_1)
				
			elif count_gapdh_grch39 > count_gapdh_grch38:
				print ('GRCH 39')
				sequencing = 'pair-end'
				count_2 = int(pysam.view("-f0X50",'-c', path, position_gapdh_grch39)) # Conversion -f80 to hexa
				count_1 = int(pysam.view("-f0X60",'-c', path, position_gapdh_grch39)) # Conversion -f96 to hexa
				
			if count_gapdh_grch38 == 0 and count_gapdh_grch39 == 0:
				sequencing = 'single-end'
				count_2_gapdh_grch38 = int(pysam.view("-f0X10",'-c', path, position_gapdh_grch38))
				count_1_gapdh_grch38 = int(pysam.view("-F0X10",'-c', path, position_gapdh_grch38))
				count_2_gapdh_grch39 = int(pysam.view("-f0X10",'-c', path, position_gapdh_grch39))
				count_1_gapdh_grch39 = int(pysam.view("-F0X10",'-c', path, position_gapdh_grch39))
				count_1 = count_1_gapdh_grch39
				count_2 = count_2_gapdh_grch39
				if (count_2_gapdh_grch38 + count_1_gapdh_grch38) > (count_2_gapdh_grch39+count_1_gapdh_grch39):
					count_1 = count_1_gapdh_grch38
					count_2 = count_2_gapdh_grch38
		else:
			position_gapdh_grch39 = 'chr12:6,537,097-6,537,227'
			try:
				count_1 = int(pysam.view("-f1",'-c', path, position_gapdh_grch39))

				if count_1 == 0 :
					sequencing = 'single-end'
					count_1 = int(pysam.view("-f0X10",'-c', path, position_gapdh_grch39))
					count_2 = int(pysam.view("-F0X10",'-c', path, position_gapdh_grch39))
				else:
					sequencing = 'pair-end'
					count_1 = int(pysam.view("-f0X50",'-c', path, position_gapdh_grch39)) # Conversion -f80 to hexa
					count_2 = int(pysam.view("-f0X60",'-c', path, position_gapdh_grch39)) # Conversion -f96 to hexa
			except pysam.utils.SamtoolsError: 

				return '',''

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
				print ('Guessing library for this Bam file %s fail. Adding unstranded library ! ' , path)
				type_library = 'unstranded'

		return type_library, sequencing

			



