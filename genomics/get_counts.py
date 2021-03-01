import os, logging, time, subprocess, pickle, multiprocessing, os, _thread, csv, collections
import genomics.get_alignments as get_alig
import pandas as pd
from pathos.multiprocessing import ProcessPool
import utils.useful_functions as uf


NUM_WORKERS =  multiprocessing.cpu_count()
__author__ = "Maria Virginia Ruiz Cuevas"


class GetCounts:

	def __init__(self, path_to_output_folder, name_exp, mode):
		self.path_to_output_folder = path_to_output_folder+'res/'
		self.name_exp = name_exp
		self.mode = mode
		self.path_to_output_folder_alignments = path_to_output_folder+'alignments/'


	def ribo_counts(self, perfect_alignments, bam_files_list):

		df_counts = pd.DataFrame({'A' : []})
		exist = os.path.exists(self.path_to_output_folder+self.name_exp+'_ribo_count.csv')
		
		if not exist:
			t_0 = time.time()
			logging.info('Alignments on ribosome profiling information ')

			to_write = 'Peptide\tPosition\tStrand\tName_Sample\tCount\n'
			
			perfect_alignments_to_return = {} 

			keys = []
			values = []
			data = []
			info_bams = []
			bams = []

			for name_sample, info_bam in bam_files_list.items():
				info_bams.append((name_sample,info_bam))
			
			keys = perfect_alignments.keys()
			values = perfect_alignments.values()
			pool = ProcessPool(nodes=NUM_WORKERS)

			bams.extend([info_bams]*len(keys))
			results = pool.map(self.get_counts_sample, bams, keys, values)
			for res in results:
				if len(res[1]) > 0:
					to_write += res[0]
					data.extend(res[1])
					for count_align in res[1]: 
						if count_align[-1] > 0:
							key = count_align[0]+'_'+count_align[1]
							count = count_align[3]
							perfect_alignments[key][-1] = count 
							perfect_alignments_to_return[key] = perfect_alignments[key]


			if len(data) > 0 :
				df_counts = pd.DataFrame(data, columns=['Peptides', 'Alignments', 'BAM Files', 'Read Counts'])
				df_counts = df_counts.groupby(['Peptides', 'BAM Files'])['Read Counts'].sum().reset_index()
				df_counts = df_counts.pivot("Peptides", "BAM Files", "Read Counts")
				_thread.start_new_thread(self.save_info_counts, (df_counts, to_write, '_ribo_count.csv'))
			
			name_path = self.path_to_output_folder_alignments+'Alignments_ribo_information.dic'
			with open(name_path, 'wb') as handle:
				pickle.dump(perfect_alignments_to_return, handle, protocol=pickle.HIGHEST_PROTOCOL)

			name_path = self.path_to_output_folder_alignments+'Alignments_information.dic'
			with open(name_path, 'wb') as handle:
				pickle.dump(perfect_alignments, handle, protocol=pickle.HIGHEST_PROTOCOL)

			t_2 = time.time()
			total = t_2-t_0
			logging.info('Total time run function get_counts to end : %f min', (total/60.0))
		else:
			logging.info('Count information already collected in the output folder : %s --> Skipping this step!', self.path_to_output_folder+self.name_exp+'_ribo_count.csv')
			
			df_counts = pd.read_csv(self.path_to_output_folder+self.name_exp+'_ribo_count.csv', index_col=0)

			with open(self.path_to_output_folder_alignments+'Alignments_information.dic', 'rb') as fp:
				perfect_alignments = pickle.load(fp)

		return perfect_alignments, df_counts


	def get_counts(self, perfect_alignments, bam_files_list):

		df_counts = pd.DataFrame()
		df_counts_filtered = pd.DataFrame()
		exist_rna = os.path.exists(self.path_to_output_folder+self.name_exp+'_rna_count.csv')
		exist_ribo = True

		if self.mode == 'filter':
			with open(self.path_to_output_folder_alignments+'Alignments_ribo_information.dic', 'rb') as fp:
				perfect_alignments_ribo = pickle.load(fp)
			data_filtered = []
			exist_ribo = os.path.exists(self.path_to_output_folder+self.name_exp+'_rna_ribo_count.csv')

		if not exist_rna and not exist_ribo:
			t_0 = time.time()
			to_write = 'Peptide\tPosition\tStrand\tName_Sample\tCount\n'
			
			keys = []
			values = []

			#for k in sorted(perfect_alignments, key=perfect_alignments.get, reverse=True):
			#	keys.append(k)
			#	values.append(perfect_alignments[k])

			data = []
			
			info_bams = []
			bams = []
			for name_sample, info_bam in bam_files_list.items():
				info_bams.append((name_sample,info_bam))
			
			keys = perfect_alignments.keys()
			values = perfect_alignments.values()
			pool = ProcessPool(nodes=NUM_WORKERS)

			bams.extend([info_bams]*len(keys))
			results = pool.map(self.get_counts_sample, bams, keys, values)
			for res in results:
				if len(res[1]) > 0:
					to_write += res[0]
					data.extend(res[1])
					for count_align in res[1]: 
						key = count_align[0]+'_'+count_align[1]
						count = count_align[3]
						perfect_alignments[key][-2] = count
						if self.mode == 'filter':
							try:
								count_ribo = perfect_alignments_ribo[key][-1]
								if count_ribo > 0:
									data_filtered.extend(res[1])
							except KeyError:
								pass

			# Replace dic
			name_path = self.path_to_output_folder_alignments+'/Alignments_information.dic'
			with open(name_path, 'wb') as handle:
				pickle.dump(perfect_alignments, handle, protocol=pickle.HIGHEST_PROTOCOL)


			if len(data) > 0 :
				df_counts = pd.DataFrame(data, columns=['Peptides', 'Alignments', 'BAM Files', 'Read Counts'])
				df_counts = df_counts.groupby(['Peptides', 'BAM Files'])['Read Counts'].sum().reset_index()
				df_counts = df_counts.pivot("Peptides", "BAM Files", "Read Counts")
				_thread.start_new_thread(self.save_info_counts, (df_counts, to_write, '_rna_count.csv'))
				if self.mode == 'filter':
					df_counts_filtered = pd.DataFrame(data_filtered, columns=['Peptides', 'Alignments', 'BAM Files', 'Read Counts'])
					df_counts_filtered = df_counts_filtered.groupby(['Peptides', 'BAM Files'])['Read Counts'].sum().reset_index()
					df_counts_filtered = df_counts_filtered.pivot("Peptides", "BAM Files", "Read Counts")
					_thread.start_new_thread(self.save_info_counts, (df_counts_filtered, to_write, '_rna_ribo_count.csv'))
			
			t_2 = time.time()
			total = t_2-t_0
			logging.info('Total time run function get_counts to end : %f min', (total/60.0))
		else:
			logging.info('Count information already collected in the output folder : %s --> Skipping this step!', self.path_to_output_folder+self.name_exp+'_rna_count.csv')
			df_counts = pd.read_csv(self.path_to_output_folder+self.name_exp+'_rna_count.csv', index_col=0)
			
			with open(self.path_to_output_folder_alignments+'/Alignments_information.dic', 'rb') as fp:
				perfect_alignments = pickle.load(fp)

			df_counts_filtered = pd.read_csv(self.path_to_output_folder+self.name_exp+'_rna_ribo_count.csv', index_col=0)
			
		return df_counts, perfect_alignments, df_counts_filtered

	def get_counts_sample(self, bams, peptide_alignment, info_alignment):
		to_return = []
		to_print = ''
		total_for_peptide_in_bam_files = 0

		for info_bam in bams:
			count = 0
			name_sample = info_bam[0]
			bam_file = info_bam[1][0]
			library = info_bam[1][1]
			sens = info_bam[1][2]

			peptide = peptide_alignment.split('_')[0]
			alignment = peptide_alignment.split('_')[1]
			
			strand = info_alignment[0]
			sequence = info_alignment[1]
			
			chr = alignment.split(':')[0]
			region_to_query = chr+':'+alignment.split(':')[1].split('-')[0]+'-'+alignment.split(':')[1].split('-')[-1]
			
			to_print = ''
			count = self.get_depth_with_view(region_to_query, bam_file, library, sens, strand, sequence)
			total_for_peptide_in_bam_files += count

			to_print += peptide+'\t'+alignment+'\t'+strand+'\t'+name_sample+'\t'+str(count)+'\n'
			to_return.append([peptide, alignment, name_sample, count])
		 
		if total_for_peptide_in_bam_files == 0:
			return '', []
		else:
			return to_print, to_return

	def save_info_counts(self, df, to_write, type_save):
		file_to_open = open(self.path_to_output_folder+self.name_exp+type_save.split('.')[0]+'.info', 'w')
		file_to_open.write(to_write)
		file_to_open.close()
		df.to_csv(self.path_to_output_folder+self.name_exp+type_save, index=True, header=True)
		logging.info('Counts Information saved to : %s ', self.path_to_output_folder+self.name_exp+type_save)


	def get_depth_with_view(self, region_to_query, bam_file, library, sens, strand, seq):
		contReads = 0
		library = library.lower()
		sens = sens.lower()
		
		if strand == '-':
			seq = uf.reverseComplement(seq)
	
		command = self.set_command_samtools(library, sens, strand, region_to_query, bam_file, seq)
		p_1 = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
		out, err = p_1.communicate()
		contReads = int(out.strip())
		return contReads

	def set_command_samtools(self, library, sens, strand, region_to_query, bam_file, seq):

		command = ''
		if library == 'single-end':
			if ((strand == '+' and sens == 'forward') or(strand == '-' and sens == 'reverse')):
				command='samtools view -F0X10 '+bam_file+' '+region_to_query+' | grep '+seq+' | sort -k1 | uniq -f1 | wc -l'
			elif ((strand == '-' and sens == 'forward') or(strand == '+' and sens == 'reverse')):
				command='samtools view -f0X10 '+bam_file+' '+region_to_query+' | grep '+seq+' | sort -k1 | uniq -f1 | wc -l'
		else:
			if ((strand == '+' and sens == 'forward') or (strand == '-' and sens == 'reverse')):
				command = 'command_1=$(samtools view -f0X60 '+bam_file+' '+region_to_query+ ' | grep '+seq+' | sort -k1 | uniq -f1); c=$(echo "${command_1}"); command_2=$(samtools view -f0X90 '+bam_file+' '+region_to_query+ ' | grep '+seq+' | sort -k1 | uniq -f1); d=$(echo "${command_2}"); e="$c""$d" ; echo "$e" | grep '+seq +' | sort -k1 | uniq -f1 | wc -l'
			elif ((strand == '-' and sens == 'forward') or(strand == '+' and sens == 'reverse')):
				command = 'command_1=$(samtools view -f0X50 '+bam_file+' '+region_to_query+ ' | grep '+seq+' | sort -k1 | uniq -f1); c=$(echo "${command_1}"); command_2=$(samtools view -f0XA0 '+bam_file+' '+region_to_query+ ' | grep '+seq+' | sort -k1 | uniq -f1); d=$(echo "${command_2}"); e="$c""$d" ; echo "$e" | grep '+seq +' | sort -k1 | uniq -f1 | wc -l'
		
		return command


