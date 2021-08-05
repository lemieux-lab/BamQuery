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


# Modify this line in order to add the optional arguments for you Query
# For instance, if you are adding the strandedness optional argument, the command variable should look like this:
# command="bamquery $input_folder $id_query $strandedness"

command="bamquery $input_folder $id_query"


echo $command | qsub -V -l nodes=1:ppn=32,mem=50gb,walltime=48:00:00 -j oe -N $id_query -d $input_folder 


echo 'Query launched for '$id_query' Done!'




