import os, logging, time, pickle, sys, getopt, os, pysam
from os import listdir
from os.path import isfile, join
from pathos.multiprocessing import ProcessPool

__author__ = "Maria Virginia Ruiz"
__email__ = "maria.virginia.ruiz.cuevas@umontreal.ca"

path_to_lib = '/'.join(os.path.abspath(__file__).split('/')[:-3])+'/lib/'

files_with_not_permission = []


class GetPrimaryReadCountBamFiles:

	def __init__(self):
		pass

	def get_all_counts_bam_file(self, bam_file_info):
		
		path = bam_file_info
		
		name_bam_file = "_".join(path.split('/')[:-1][-2:])
		contReads = int(pysam.view("-F0X100",'-c', path).strip()) 
		
		bam_file_info_to_return = [name_bam_file, contReads]
		
		logging.info('BAM file : %s Path: %s Cont Reads: %s ', name_bam_file, path, str(contReads))
		
		return bam_file_info_to_return

	def remove_lock_to_bam_files_info_dic(self, path_to_lock_file, path_to_input_folder):
		with open(path_to_lock_file, 'r') as file:
			first_line = file.readline().strip()
		if first_line == path_to_input_folder:
			first_line = ''
			with open(path_to_lock_file, 'w') as file:
				file.writelines(first_line) 
			file.close()
		else:
			print ('Error in the path_to_lock_file: contents ', first_line, ' and the ', path_to_input_folder, ' not the same.')
		

	def grant_access_to_bam_files_info_dic(self, path_to_lock_file, path_to_input_folder, timeout=300):
		start_time = time.time()
		with open(path_to_lock_file, 'r') as file:
			first_line = file.readline().strip()
		while (time.time() - start_time < timeout) and first_line != path_to_input_folder:
			try:
				with open(path_to_lock_file, 'r') as file:
					first_line = file.readline().strip()
				if first_line == '':
					with open(path_to_lock_file, 'w') as file:
						file.write(path_to_input_folder)
					file.close()
					return True
				time.sleep(1)
			except IOError:
				pass
		raise Exception ("\nTimeout to grant access to Bam_files_info.dic. Another BamQuery process may be using the dictionary." )


	def set_values(self, path_to_lib, bam_files_list, threads):

		path_to_all_counts_file = path_to_lib+"Bam_files_info.dic"
		path_to_lock_file = path_to_lib+"lock_dic"

		logging.info('Total Bam Files to Query : %d %s', len(bam_files_list), bam_files_list)
		
		if len(bam_files_list) >= threads:
			pool = ProcessPool(nodes=threads)
		else:
			pool = ProcessPool(nodes=len(bam_files_list))
		
		results = pool.map(self.get_all_counts_bam_file, bam_files_list)

		if self.grant_access_to_bam_files_info_dic(path_to_lock_file, 'getting_primary_read_count_bam_files'):
			logging.info('Lock access to Bam_files_info.dic')
			pass
		else:
			raise Exception("\nCould not acquire the lock to consult the Bam_files_info.dic. Check if there is another run of BamQuery that wrote to the lock file and it is still running." )

		with open(path_to_all_counts_file, 'rb') as fp:
			dictionary_total_reads_bam_files = pickle.load(fp)

		for bam_file_info_to_return in results:
			name_bam_file = bam_file_info_to_return[0]
			count = bam_file_info_to_return[1]
			dictionary_total_reads_bam_files[name_bam_file][1] = int(count)

		with open(path_to_all_counts_file, 'wb') as handle:
			pickle.dump(dictionary_total_reads_bam_files, handle, protocol=pickle.HIGHEST_PROTOCOL)

		logging.info('Primary read counts for %s have been saved in the Bam_files_info.dic.', str(bam_files_list))

		self.remove_lock_to_bam_files_info_dic(path_to_lock_file, 'getting_primary_read_count_bam_files')
		logging.info('Unlock Bam_files_info')

		
def main(argv):

	path_to_input_folder = ''
	path_to_output_folder = ''

	try:
		opts, args = getopt.getopt(argv,"hi:o:t:",["path_to_input_folder=", "path_to_output_folder=", "threads="])
	except getopt.GetoptError:
		print ('GetPrimaryReadCountBamFiles.py -i <path_to_input_folder> -o <path_to_output_folder> -t <threads>')
		sys.exit(2)

	for opt, arg in opts:
		if opt == '-h':
			print ('GetPrimaryReadCountBamFiles.py -i <path_to_input_folder>')
			sys.exit()
		elif opt in ("-i", "--path_to_input_folder"):
			path_to_input_folder = arg
		elif opt in ("-o", "--path_to_output_folder"):
			path_to_output_folder = arg
		elif opt in ("-t", "--threads"):
			threads = int(arg)

	
	global path_to_lib
	path_to_lib = '/'.join(os.path.abspath(__file__).split('/')[:-3])+'/lib/'
	
	nameLog = path_to_output_folder+'logs/Get_Read_Count_BAM_directories.log'
	logging.basicConfig(filename=nameLog, filemode='w', level=logging.INFO, format='%(asctime)s %(message)s')

	t0 = time.time()
	Get_Read_Count_obj = GetPrimaryReadCountBamFiles()
	
	with open(path_to_input_folder, 'rb') as fp:
		bam_files_list = pickle.load(fp)
	
	Get_Read_Count_obj.set_values(path_to_lib, bam_files_list, threads)
	t2 = time.time()
	total = t2-t0
	logging.info('Total time run Primary read count to end : %f min', (total/60.0))

if __name__ == "__main__":
	main(sys.argv[1:])

