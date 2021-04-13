import os, logging, threading, time, subprocess, concurrent.futures, getpass, pickle, sys, getopt, os, pysam
from os import listdir
from os.path import isfile, join

__author__ = "Maria Virginia Ruiz"
__email__ = "maria.virginia.ruiz.cuevas@umontreal.ca"


class GetPrimaryReadCountBamFiles:

	def __init__(self):
		pass

	def get_all_counts(self, path_to_input_folder, path_to_output_folder, mode, strandedness):

		path_to_lib = '/'.join(os.path.abspath(__file__).split('/')[:-3])+'/lib/'
		path_to_all_counts_file = path_to_lib+"allcounts.dic"
		exist = os.path.exists(path_to_all_counts_file)
		
		if exist :
			with open(path_to_all_counts_file, 'rb') as fp:
				dictionary_total_reads_bam_files = pickle.load(fp)
		else:
			dictionary_total_reads_bam_files = {}

		self.bam_files_list = {}
		self.bam_ribo_files_list = {}

		try:
			bam_files = path_to_input_folder+'BAM_directories.tsv'
			self.bam_files_list = self.get_info_bamfiles(bam_files, strandedness, path_to_output_folder)

		except FileNotFoundError:
			logging.info('Bam File %s doesn\'t exist ', path_to_input_folder+'BAM_directories.tsv')

		if mode == 'translation':
			try:
				bam_files = path_to_input_folder+'BAM_Ribo_directories.tsv'
				self.bam_ribo_files_list = self.get_info_bamfiles(bam_files, strandedness, path_to_output_folder)

			except FileNotFoundError:
				logging.info('If running translation mode you must include a list of Ribo Bam Files. Bam File %s doesn\'t exist ', path_to_input_folder+'BAM_Ribo_directories.tsv')


	def get_info_bamfiles(self, bam_files, strandedness, path_to_output_folder):

		bam_files_list = {}

		with open(bam_files) as f:

			for index, line in enumerate(f):
				line = line.strip().split('\t')
				name = line[0]
				path = line[1]
				
				if '.bam' not in path:
					bam_file_paths = self.search_bam_files(path)

					for bam_file in bam_file_paths:
						if strandedness:
							library, sequencing = self.get_type_library(bam_file)
						else:
							library = 'unstranded' 
							sequencing = 'unstranded' 
						name_bam_file = "_".join(bam_file.split('/')[:-1][-2:])
						bam_files_list[name_bam_file] = [bam_file, sequencing, library, name]
				else:
					if strandedness:
						library, sequencing = self.get_type_library(bam_file)
					else:
						library = 'unstranded' 
						sequencing = 'unstranded' 
					name_bam_file = "_".join(path.split('/')[:-1][-2:])
					bam_files_list[name_bam_file] = [path, sequencing, library, name]
		
		get_read_counts_path = os.path.abspath(__file__)
		command = 'python '+get_read_counts_path+' -i '+bam_files+' -o '+ path_to_output_folder
		subprocess.Popen(command, shell=True, stdout=subprocess.DEVNULL, close_fds=True)
		logging.info('Total Bam Files to Query : %d', len(bam_files_list))
		return bam_files_list

	def get_all_counts_bam_file(self, bam_file_info):
		
		path = bam_file_info
		name_bam_file = "_".join(path.split('/')[:-1][-2:])
		contReads = -1
		count = pysam.view("-F0X100",'-c', name_bam_file) 
		bam_file_info_to_return = [name_bam_file, path, contReads]
		logging.info('BAM file : %s Path: %s Cont Reads: %s ', name_bam_file, path, str(contReads))

		return bam_file_info_to_return


	def set_values(self, path_to_lib, bam_files_dic):

		path_to_all_counts_file = path_to_lib+"allcounts.dic"
		exist = os.path.exists(path_to_all_counts_file)
		
		if exist :
			with open(path_to_all_counts_file, 'rb') as fp:
				dictionary_total_reads_bam_files = pickle.load(fp)
		else:
			dictionary_total_reads_bam_files = {}

		list_bams_to_assess = []

		some_change = False

		for name_bam_file, bam_file_info in bam_files_dic.items():
			path = bam_file_info
			
			try:
				info_bam_file = dictionary_total_reads_bam_files[name_bam_file]
				path_saved = info_bam_file[0]
				count_reads_saved = info_bam_file[1]
				
				if path_saved != path:
					
					if  path_saved == '':
						some_change = True
						info_bam_file[0] = path
					else:
						logging.info('BAM file is already in the dictionary, however the path for the BAM file is not the same. \
							Path found: %s  path in BAM_directories.tsv or BAM_Ribo_directories.tsv is %s ', path_saved, path)
				else:
					logging.info('%s is already in the dictionary, total reads : %s ', path, str(count_reads_saved))
			except KeyError:
				list_bams_to_assess.append(path)

		user = getpass.getuser()
		logging.info('Total Bam Files to Query : %d %s', len(list_bams_to_assess), list_bams_to_assess)
		
		if some_change:
			with open(path_to_lib+"allcounts.dic", 'wb') as handle:
				pickle.dump(dictionary_total_reads_bam_files, handle, protocol=pickle.HIGHEST_PROTOCOL)

		if len(list_bams_to_assess) > 0:
			with concurrent.futures.ThreadPoolExecutor() as executor:	
				futures = [executor.submit(self.get_all_counts_bam_file, param) for param in list_bams_to_assess]
			
			for f in futures:
				bam_file_info_to_return = f.result()
				dictionary_total_reads_bam_files[bam_file_info_to_return[0]] = [bam_file_info_to_return[1], bam_file_info_to_return[2], user]

			with open(path_to_lib+"allcounts.dic", 'wb') as handle:
				pickle.dump(dictionary_total_reads_bam_files, handle, protocol=pickle.HIGHEST_PROTOCOL)
			logging.info('Saved information of read count for each BAM file')
		else:
			logging.info('All bam files into allcounts dictionary !')


	def search_bam_files(self, path):
		try:
			bam_file_paths = []
			components = [f for f in listdir(path)]
			
			for component in components:
				if component[0] != '.' and component[0] != '_':
					path_component = join(path, component)
					if isfile(path_component) :
						if path_component[-4:] == '.bam':
							bam_file_paths.append(path_component)
					else:
						bam_file_paths = self.search_bam_file(bam_file_paths, path_component)
			
			return bam_file_paths
		except PermissionError:
			logging.info('You don\'t have permission to read the files contained in %s !', path)

	def search_bam_file(self, bam_file_paths, path): 
		try:
			components = [f for f in listdir(path)]
			for component in components:
				
				if component[0] != '.' and component[0] != '_':
					path_component = join(path, component)
					
					if isfile(path_component) :
						if path_component[-4:] == '.bam':
							bam_file_paths.append(path_component)
					else:
						self.search_bam_file(bam_file_paths, path_component)
			return bam_file_paths
		except PermissionError:
			logging.info('You don\'t have permission to read the files contained in %s !', path)


	def get_type_library(self, path):

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
				logging.info('Guessing library for this Bam file %s fail. Adding forward library ! ' , path)
				type_library = 'forward'

		return type_library, sequencing


def main(argv):

	path_to_input_folder = ''
	path_to_output_folder = ''

	try:
		opts, args = getopt.getopt(argv,"hi:o:",["path_to_input_folder=", "path_to_output_folder="])
	except getopt.GetoptError:
		print ('GetPrimaryReadCountBamFiles.py -i <path_to_input_folder> -o <path_to_output_folder>')
		sys.exit(2)

	for opt, arg in opts:
		if opt == '-h':
			print ('GetPrimaryReadCountBamFiles.py -i <path_to_input_folder>')
			sys.exit()
		elif opt in ("-i", "--path_to_input_folder"):
			path_to_input_folder = arg
		elif opt in ("-o", "--path_to_output_folder"):
			path_to_output_folder = arg

	
	global path_to_lib
	path_to_lib = '/'.join(os.path.abspath(__file__).split('/')[:-3])+'/lib/'
	
	name_log_file = path_to_input_folder.split('/')[-1].split('.')[0]
	nameLog = path_to_output_folder+'logs/Get_Read_Count_'+name_log_file+'.log'
	logging.basicConfig(filename=nameLog, filemode='w', level=logging.INFO, format='%(asctime)s %(message)s')

	t0 = time.time()
	bam_files_dic = {}
	bam_files = path_to_input_folder
	Get_Read_Count_obj = GetPrimaryReadCountBamFiles()

	with open(bam_files) as f:

		for index, line in enumerate(f):
			line = line.strip().split('\t')
			path = line[1]
			
			if '.bam' not in path:
				bam_file_paths = Get_Read_Count_obj.search_bam_files(path)

				for bam_file in bam_file_paths:
					name_bam_file = "_".join(bam_file.split('/')[:-1][-2:])
					bam_files_dic[name_bam_file] = bam_file
			else:
				name_bam_file = "_".join(path.split('/')[:-1][-2:])
				bam_files_dic[name_bam_file] = path

	Get_Read_Count_obj.set_values(path_to_lib, bam_files_dic)
	t2 = time.time()
	total = t2-t0
	logging.info('Total time run function BamQuery to end : %f min', (total/60.0))

if __name__ == "__main__":
	main(sys.argv[1:])

