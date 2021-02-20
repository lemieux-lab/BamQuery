import os, logging, threading, time, subprocess, concurrent.futures, getpass, pickle, sys, getopt

__author__ = "Maria Virginia Ruiz"

class GetPrimaryReadCountBamFiles:

	def __init__(self):
		pass


	def get_all_counts(self, path_to_input_folder, path_to_output_folder):

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

			with open(bam_files) as f:
				for index, line in enumerate(f):
					line = line.strip().split('\t')
					name_bam_file = "_".join(line[1].split('/')[:-1][-2:])
					self.bam_files_list[name_bam_file] = [line[1], line[2], line[3], line[0]]
			
			get_read_counts_path = os.path.abspath(__file__)
			command = 'python '+get_read_counts_path+' -i '+bam_files+' -o '+ path_to_output_folder
			subprocess.Popen(command, shell=True, stdout=subprocess.DEVNULL, close_fds=True)
			logging.info('Total Bam Files to Query : %d', len(self.bam_files_list))

		except FileNotFoundError:
			logging.info('Bam File %s doesn\'t exist ', path_to_input_folder+'BAM_directories.tsv')

		try:
			bam_files = path_to_input_folder+'BAM_Ribo_directories.tsv'

			with open(bam_files) as f:
				for index, line in enumerate(f):
					line = line.strip().split('\t')
					name_bam_file = "_".join(line[1].split('/')[:-1][-2:])
					self.bam_ribo_files_list[name_bam_file] = [line[1], line[2], line[3], line[0]]
			
			get_read_counts_path = os.path.abspath(__file__)
			command = 'python '+get_read_counts_path+' -i '+bam_files+' -o '+ path_to_output_folder
			subprocess.Popen(command, shell=True, stdout=subprocess.DEVNULL, close_fds=True)
			logging.info('Total Ribo Bam Files to Query : %d', len(self.bam_ribo_files_list))

		except FileNotFoundError:
			logging.info('If running filter mode you must include a list of Ribo Bam Files. ')


	def get_all_counts_bam_file(self, bam_file_info):
		
		path = bam_file_info[0]
		name_bam_file = "_".join(path.split('/')[:-1][-2:])
		contReads = -1
		command = 'samtools view -F 256 '+path+' | wc -l'
		p_1 = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
		out, err = p_1.communicate()
		contReads = int(out.strip())
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

		for name_bam_file, bam_file_info in bam_files_dic.items():
			path = bam_file_info[0]
			
			try:
				info_bam_file = dictionary_total_reads_bam_files[name_bam_file]
				path_saved = info_bam_file[0]
				count_reads_saved = info_bam_file[1]
				if info_bam_file[0] != path: 
					logging.info('BAM file is already in the DB, however the path for the BAM file is not the same. \
						Path found: %s  path in BAM_directories.tsv %s ', path_saved, path)
				else:
					logging.info('%s is already in the DB, total reads : %s ', path, str(count_reads_saved))
			except KeyError:
				list_bams_to_assess.append(bam_file_info)

		user = getpass.getuser()
		logging.info('Total Bam Files to Query : %d %s', len(list_bams_to_assess), list_bams_to_assess)

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
	
	nameLog = path_to_output_folder+'logs/Get_Read_Count_Res.log'
	logging.basicConfig(filename=nameLog, filemode='w', level=logging.INFO, format='%(asctime)s %(message)s')

	t0 = time.time()
	bam_files_dic = {}
	bam_files = path_to_input_folder
	
	with open(bam_files) as f:
		for index, line in enumerate(f):
			line = line.strip().split('\t')
			name_bam_file = "_".join(line[1].split('/')[:-1][-2:])
			bam_files_dic[name_bam_file] = [line[1], line[2], line[3], line[0]]

	Get_Read_Count_obj = GetPrimaryReadCountBamFiles()

	Get_Read_Count_obj.set_values(path_to_lib, bam_files_dic)
	t2 = time.time()
	total = t2-t0
	logging.info('Total time run function BamQuery to end : %f min', (total/60.0))

if __name__ == "__main__":
	main(sys.argv[1:])

