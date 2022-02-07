#!/bin/bash


export REF_PATH=~/tmp
module load R/4.0.0
module load bamquery
module add torque

echo '#################### Launching BamQuery ######################'


# positional arguments:
input_folder=''							# Path to the input folder where to find BAM_directories.tsv and peptides.tsv
id_query=''								# BamQuery search Id

# optional arguments:
mode=''									# --mode MODE           BamQuery search mode : normal / translation. Normal is the mode by default.
										# The variable mode should look like this if you specified it : mode='--mode normal' 

strandedness=''							# --strandedness        Take into account strandedness of the samples
										# The variable strandedness should look like this if you specified it : strandedness='--strandedness' 

th_out=''								# --th_out th_out 		Threshold to assess expression comparation with other tissues
										# The variable th_out should look like this if you specified it : th_out='--th_out 8.55'

light=''								# --light               Display only the count and norm count for peptides and regions
										# The variable light should look like this if you specified it : light='--light'

dbSNP=''								# The dbSNP variable allows you to choose between three versions of dbSNPs: 149 / 151 / 155. dbSNP 149 is the default.
										# If you don't want to use any release specify 0 for this argument. Example: dbSNP='--dbSNP 151'

c=''									# This option allows to choose between the most COMMON SNPs from the dbSNP release that you choose with the argument above.
										# The variable mode should look like this if you specified it : c='--c'

plots=''								# This option sets BamQuery to produce pie charts in the biotype analysis step. 
										# The variable mode should look like this if you specified it : plots='--plots'

genome_version=''						# This option allows to choose between three genome versions : v26_88 / v33_99 / v38_104. genome version v26_88 is the default.
										# The variable mode should look like this if you specified it : genome_version='--genome_version v26_88'

maxmm=''								# This option causes the STAR aligner to find all possible locations in the genome for MCSs. It is advisable to use this option for ERE peptides.
										# The variable should look like this: maxmm='--maxmm'

# Modify this line in order to add the optional arguments for you Query
# For instance, if you are adding the strandedness optional argument, the command variable should look like this:
# command="bamquery $input_folder $id_query $strandedness"

command="bamquery $input_folder $id_query"


echo $command | qsub -V -l nodes=1:ppn=32,mem=50gb,walltime=48:00:00 -j oe -N $id_query -d $input_folder 


echo 'Query launched for '$id_query' Done!'




