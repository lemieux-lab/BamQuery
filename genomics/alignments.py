import os, logging, time, subprocess, pickle, multiprocessing, os, _thread, csv
import genomics.get_alignments as get_alig
import collections

__author__ = "Maria Virginia Ruiz Cuevas"
__email__ = "maria.virginia.ruiz.cuevas@umontreal.ca"

NUM_WORKERS =  multiprocessing.cpu_count()

class Alignments:

	def __init__(self, path_to_output_folder, name_exp):
		self.path_to_output_folder_genome_alignments = path_to_output_folder+'genome_alignments/'
		self.path_to_output_folder_alignments = path_to_output_folder+'alignments/'
		self.name_exp = name_exp
		self.path_to_lib = '/'.join(os.path.abspath(__file__).split('/')[:-3])+'/lib/'

	def alignment_cs_to_genome(self, set_peptides):

		exist = os.path.exists(self.path_to_output_folder_genome_alignments+'/Aligned.out.sam')

		if not exist:
			t_0 = time.time()
			inputFilesR1_1 = self.path_to_output_folder_genome_alignments+self.name_exp+'.fastq'
			genomeDirectory = self.path_to_lib+'/Index_BAM_Query/'
			outputFilterMatchInt = 20
			seed = 20
			anchor = 2000										
			maxMulti = 2000 						
			dbSNPFile = self.path_to_lib+'dbSNP149_all.vcf'

			command='module add star/2.7.1a; STAR --runThreadN '+ str(NUM_WORKERS)+\
					' --genomeDir '+genomeDirectory+' --seedSearchStartLmax '+str(seed)+\
					' --alignEndsType EndToEnd --sjdbOverhang 32 --sjdbScore 2 --alignSJDBoverhangMin 3 --alignSJoverhangMin 20 --outFilterMismatchNmax 4 --winAnchorMultimapNmax '+\
					str(anchor)+' --outFilterMultimapNmax '+str(maxMulti)+' --outFilterMatchNmin '+str(outputFilterMatchInt)+' --genomeConsensusFile '+\
					dbSNPFile+' --readFilesIn  '+inputFilesR1_1+' --outSAMattributes NH HI MD --outFileNamePrefix '+self.path_to_output_folder_genome_alignments # 3286 #3276
			
			logging.info('Command to Align using STAR : %s ', command)
			p_1 = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
			out, err = p_1.communicate()
			t_2 = time.time()
			total = t_2-t_0
			print ("Total time run function alignment_CS_to_genome End : %s min" % (total/60.0))
			logging.info('Total time run function alignment_CS_to_genome to end : %f min', (total/60.0))
		else:
			logging.info('Alignment file already exists in the output folder : %s --> Skipping this step!', self.path_to_output_folder_alignments+'/Aligned.out.sam')
		
		perfect_alignments = self.get_alignments(set_peptides)
		
		return perfect_alignments

	def get_alignments(self, set_peptides):
		t_0 = time.time()
		sam_file = self.path_to_output_folder_genome_alignments+'/Aligned.out.sam'
		
		exist = os.path.exists(self.path_to_output_folder_alignments+'/Alignments_information.dic')
		perfect_alignments = {}

		if not exist:
			res_star = get_alig.get_alignments(sam_file)
			
			t_2 = time.time()
			total = t_2-t_0
			print ("Total time run function get_alignments End : %s " % (total/60.0))
			logging.info('Total time run function get_alignments to end : %f min', (total/60.0))
			logging.info('Total perfect aligments : %s ', str(len(res_star[0])))
			
			perfect_alignments = res_star[0]

			peptides_with_alignments = res_star[3]

			_thread.start_new_thread(self.save_output_info, (res_star[0], '_info_perfect_alignments.csv',))
			_thread.start_new_thread(self.save_output_info, (res_star[1], '_info_variants_alignments.csv',))
			_thread.start_new_thread(self.save_output_info, (res_star[2], '_info_out_alignments.csv',))
			_thread.start_new_thread(self.save_info, (res_star[0], ) )

		else:
			logging.info('Alignment information already collected in the output folder : %s --> Skipping this step!', self.path_to_output_folder_alignments+'/Alignments_information.dic')
			files_already_collected = [os.path.exists(self.path_to_output_folder_alignments+self.name_exp+'_info_perfect_alignments.csv'),
										os.path.exists(self.path_to_output_folder_alignments+self.name_exp+'_info_variants_alignments.csv'),
										os.path.exists(self.path_to_output_folder_alignments+self.name_exp+'_info_out_alignments.csv')]

			with open(self.path_to_output_folder_alignments+'/Alignments_information.dic', 'rb') as fp:
				perfect_alignments = pickle.load(fp)

			logging.info('Total perfect aligments : %s ', str(len(perfect_alignments)))
			
			peptides_with_alignments = set()
			
			try:
				with open(self.path_to_output_folder_alignments+'alignments/missed_peptides.info') as f:
					for index, line in enumerate(f):
						peptide = line.strip()
						peptides_with_alignments.add(peptide)
			except :
				peptides_keys = list(perfect_alignments.keys())
				for key in peptides_keys:
					peptides_with_alignments.add(key.split('_')[0])


			if sum(files_already_collected) > 0 :

				if not files_already_collected[0]:
					logging.info('Total perfect aligments : %s ', str(len(res_star)))
					_thread.start_new_thread(self.save_output_info, (res_star, '_info_perfect_alignments.csv',))
				# if not files_already_collected[1]:
				# 	_thread.start_new_thread(self.save_output_info, (res_star[1], '_info_variants_alignments.csv',))
				# if not files_already_collected[2]:
				# 	_thread.start_new_thread(self.save_output_info, (res_star[2], '_info_out_alignments.csv',))

		return perfect_alignments, peptides_with_alignments

	def save_info(self, res_star):
		name_path = self.path_to_output_folder_alignments+'/Alignments_information.dic'
		with open(name_path, 'wb') as handle:
			pickle.dump(res_star, handle, protocol=pickle.HIGHEST_PROTOCOL)
		logging.info('Alignments Information save to : %s ', name_path)


	def save_output_info(self, alignments_input, name_file):
		row_list = [["Peptide", "Strand", "Alignment", "MCS", "Peptide in Reference", "Diff_AA", "SNVs", "Mutations_Non_Annotated"]]

		alignments = collections.OrderedDict(sorted(alignments_input.items()))
		
		for peptide_info, info_alignment in alignments.items():
			alignment = peptide_info.split('_')[1]
			peptide = peptide_info.split('_')[0]

			strand = info_alignment[0]
			MCS = info_alignment[1]
			peptide_with_snps_local_reference = info_alignment[2]
			snvs_write = ''
			dif_aa_write = ''
			mutations_write = ''

			try:
				snvs = info_alignment[3][0]
				dif_aas = info_alignment[3][1]
				mutations_non_annotated = info_alignment[3][2]
				for snv in snvs:
					#T->G|snp:rs12036323|GenPos:239044335|MCSPos:10
					snv = '['+snv[0]+'->'+snv[1]+'|snp:'+snv[2]+'|GenPos:'+str(snv[3])+'|MCSPos:'+str(snv[4])+']'
					snvs_write += snv
				for dif_aa in dif_aas:
					dif_aa_write += '['+dif_aa+']'
				for mutation in mutations_non_annotated:
					mutation_write = '['+mutation+']'
					mutations_write += mutation_write
			except IndexError:
				pass
			row_list.append([peptide, strand, alignment, MCS, peptide_with_snps_local_reference, dif_aa_write, snvs_write, mutations_write])
		
		with open(self.path_to_output_folder_alignments+self.name_exp+name_file, 'w') as file:
			writer = csv.writer(file)
			writer.writerows(row_list)
		logging.info('Info perfect alignments saved to : %s ', self.path_to_output_folder_alignments+self.name_exp+name_file)


