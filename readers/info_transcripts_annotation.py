import os, logging, time, pickle, threading, csv, math, pysam, subprocess, sys


path_to_lib = '/'.join(os.path.abspath(__file__).split('/')[:-3])+'/lib/'

annotation = path_to_lib + 'gencode.v26.primary_assembly.annotation.gtf'

__author__ = "Maria Virginia Ruiz Cuevas"


class InfoTranscripts:

	def __init__(self):
		pass

	def set_values(self):

		exist = os.path.exists(path_to_lib+'Info_Transcripts_Annotations.dic')

		if not exist:
			get_info_transcripts_path = os.path.abspath(__file__)
			command = 'python '+get_info_transcripts_path
			subprocess.Popen(command, shell=True, stdout=subprocess.DEVNULL, close_fds=True) 


	def get_info_transcripts(self):

		self.info_transcript = {}
		transcript_count = 0
		transcript_cds = ''

		with open(annotation) as f:
			for index, line in enumerate(f):
				if '#' not in line:
					type_, chr, start, end, strand, gene_id, transcript_id, gene_type, gene_name, transcript_type, tsl = self.get_info_line(line)
					
					if type_ == 'transcript':
						try:
							infoTrans = self.info_transcript[transcript_id]
						except KeyError:
							transcript_count += 1
							# 0: chr, 1: start, 2: end, 3: strand, 4: gene_id, 5: gene_type
							# 6: gene_name, 7: transcript_id, 8: transcript_type, 9: tsl,
							# 10: start_codon, 11: stop_codon, 12: rangeSC, 13: lenght AA


							## Old
							# 0: chr, 1: start, 2: end, 3: transcript_type, 
							# 4: start_codon, 5: stop_codon, 6: gene, 
							# 7: rangeSC, 8: lenght AA, 9: tsl, 10: gene_name, 11: strand


							self.info_transcript[transcript_id] = { 'Info': [chr, start, end, strand, gene_id, gene_type, 
																			gene_name, transcript_id, transcript_type, tsl,
																			[], [], [], 0], 
																			'Introns': [], 'CDS': [],'Exons': [], 
																			'3UTR': [], '5UTR': []}
						
						# This is to get the first codon from the transcripts that don't have 
						# an start codon annotated
						if transcript_count > 1 and len(self.info_transcript[transcript_old]['CDS']) > 0 :
							cds_s = self.info_transcript[transcript_old]['CDS']
							
							# This is to get the lenght of the protein
							count_ntds = 0 
							for cds in cds_s:
								start = cds[0]
								end = cds[1]
								count_ntds += end - start + 1

							# Set length of the protein
							self.info_transcript[transcript_old]['Info'][13] = count_ntds/3
							
							if transcript_cds[-1]== '|':
								transcript_cds = self.info_transcript[transcript_old]['Info'][0]+':'+transcript_cds[:-1]
								self.info_transcript[transcript_old]['Info'].append(transcript_cds)
								transcript_cds = ''
								
							# To get the first 3 nucleotides from the transcript
							start_cds_1 = self.info_transcript[transcript_old]['CDS'][0][0]
							end_cds_1 = self.info_transcript[transcript_old]['CDS'][0][1]
							
							count = end_cds_1 - start_cds_1 + 1
							range_cds = range(start_cds_1, end_cds_1+1)
							
							if len(range_cds) < 3:
								sc = range_cds
								self.info_transcript[transcript_old]['Info'][12].append((sc[0],sc[-1]))
								start_cds_2 = self.info_transcript[transcript_old]['CDS'][1][0]
								end_cds_2 = self.info_transcript[transcript_old]['CDS'][1][1]
								range_cds = range(start_cds_2, end_cds_2+1)
								
								if strand_old == '+':
									sc = range_cds[0:3-(count)]
								else:
									sc = range_cds[-3+count:]
								self.info_transcript[transcript_old]['Info'][12].append((sc[0],sc[-1]))
							
							else:
								if strand_old == '+':
									sc = range_cds[0:3]
								else:
									sc = range_cds[-3:]
								self.info_transcript[transcript_old]['Info'][12].append((sc[0],sc[-1]))
								
					if type_ == 'exon':

						if len(self.info_transcript[transcript_id]['Exons']) > 0:
							if strand == '-':
								old_start = self.info_transcript[transcript_id]['Exons'][-1][0]-1
								self.info_transcript[transcript_id]['Introns'].append((end+1, old_start))
							else:
								old = self.info_transcript[transcript_id]['Exons'][-1][1]+1
								self.info_transcript[transcript_id]['Introns'].append((old, start-1))
						self.info_transcript[transcript_id]['Exons'].append((start,end))
						
								
					elif type_ == 'CDS':
						self.info_transcript[transcript_id]['CDS'].append((start,end))
						transcript_cds+=str(start)+'-'+str(end)+'|'

					elif  type_ == 'start_codon':
						self.info_transcript[transcript_id]['Info'][10].append((start,end))
					
					elif  type_ == 'stop_codon':
						self.info_transcript[transcript_id]['Info'][11].append((start,end))
					
					elif  type_ == 'UTR':
						
						if strand == '-':
							if start > self.info_transcript[transcript_id]['CDS'][0][1]:
								self.info_transcript[transcript_id]['5UTR'].append((start,end))
							else:
								self.info_transcript[transcript_id]['3UTR'].append((start,end))
						if strand == '+':
							if end < self.info_transcript[transcript_id]['CDS'][0][0]:
								self.info_transcript[transcript_id]['5UTR'].append((start,end))
							else:
								self.info_transcript[transcript_id]['3UTR'].append((start,end))

					if type_ in ['UTR', 'stop_codon', 'start_codon', 'CDS', 'exon', 'transcript']:
						strand_old = strand 
						transcript_old = transcript_id
						
		with open(path_to_lib+'Info_Transcripts_Annotations.dic', 'wb') as handle:
			pickle.dump(self.info_transcript, handle, protocol=pickle.HIGHEST_PROTOCOL)


	def get_info_line(self, line):

		#chrX	HAVANA	transcript	155466540	155494110	.	+	.	gene_id "ENSG00000224533.4"; transcript_id "ENST00000433624.1"; gene_type "antisense"; gene_name "TMLHE-AS1"; transcript_type "antisense"; transcript_name "TMLHE-AS1-002"; level 2; transcript_support_level "1"; tag "basic"; havana_gene "OTTHUMG00000022672.1"; havana_transcript "OTTHUMT00000058815.1";	26
		splitLine = line.strip().split("\t")
		chr = splitLine[0]
		type_ = splitLine[2]
		start = splitLine[3]
		end = splitLine[4]
		strand = splitLine[6]
		gene_id = splitLine[8].split('gene_id')[1].strip().split("\"")[1]
		try:
			transcript_id = splitLine[8].split('transcript_id')[1].strip().split("\"")[1]
		except IndexError:
			transcript_id = splitLine[8].split(';')[1].split("\"")[1]

		gene_type = ''
		gene_name = ''

		if 'gene_type' in line:
			gene_type = splitLine[8].split('gene_type')[1].strip().split("\"")[1]

		if 'gene_name' in line:
			gene_name = splitLine[8].split('gene_name')[1].strip().split("\"")[1]

		try:
			transcript_type = line.split('transcript_type')[1].strip().split("\"")[1]
		except IndexError:
			transcript_type = gene_type

		if type_ == 'transcript':
			try:
				tsl = splitLine[8].split('transcript_support_level')[1].strip().split("\"")[1]
			except IndexError:
				tsl = 'Not defined'
		else:
			tsl = 'exon'
		return  type_, chr, int(start), int(end), strand, gene_id, transcript_id, gene_type, gene_name, transcript_type, tsl


def main(argv):

	info_transcript = InfoTranscripts()
	info_transcript.get_info_transcripts()
	

if __name__ == "__main__":
	main(sys.argv[1:])


