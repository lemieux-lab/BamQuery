====================
Command Line Options
====================

At the command line::

    BamQuery --help


.. code::

	    usage: BamQuery.py [-h] [--mode MODE] [--genome_version GENOME_VERSION]
                   [--th_out TH_OUT] [--dbSNP DBSNP] [--strandedness]
                   [--light] [--c] [--sc] [--var] [--maxmm] [--overlap]
                   [--plots] [--dev]
                   path_to_input_folder name_exp

		======== BamQuery ========

		positional arguments:
		  path_to_input_folder  Path to the input folder where to find
		                        BAM_directories.tsv and peptides.tsv
		  name_exp              BamQuery search Id

		optional arguments:
		  -h, --help            show this help message and exit
		  --mode MODE           BamQuery search mode : normal / translation
		  --genome_version GENOME_VERSION
		                        Genome version supported : v26_88 / v33_99 / v38_104
		  --th_out TH_OUT       Threshold to assess expression comparation with other
		                        tissues
		  --dbSNP DBSNP         BamQuery dbSNP : 149 / 151 / 155 / 0
		  --strandedness        Take into account strandedness of the samples
		  --light               Display only the count and norm count for peptides and
		                        regions
		  --c                   Take into account the only common SNPs from the dbSNP
		                        database chosen
		  --sc                  Query Single Cell Bam Files
		  --var                 Keep Variants Alignments
		  --maxmm               Keep High Amount Alignments
		  --overlap             Count overlapping reads
		  --plots               Plot biotype pie-charts
		  --dev                 Save all temps files

====================



Positional Arguments
====================

Let's take a look to the positional arguments:

path_to_input_folder
--------------------

To run BamQuery, you must create a folder called Input. In this folder you must include two files: BAM_directories.tsv and peptides.tsv.

.. code::

	
	Input
	----|
	    |---:> BAM_directories.tsv
	    |---:> peptides.tsv
	    

The BAM_directories.tsv file must contain the BAM/CRAM files you want to query.
The peptides.tsv file must contain the peptides for which you want to compute their expression in the BAM/CRAM files.

name_exp
--------

To identify your BamQuery run, you must include a name.


Optional Arguments
==================


*--mode*
--------

BamQuery has to modes of search : normal / translation

If no --mode argument is specified, BamQuery will run by default in the normal mode. 


	Normal mode:
	

	In normal mode, BamQuery will expect to find in the input folder `path_to_input_folder` the files BAM_directories.tsv and peptides.tsv. In this mode, BamQuery will look for the peptide locations in the BAM/CRAM file(s) in the BAM directories. Along with the expression of the heatmaps of each peptide for each BAM/CRAM file/files, you will find the biotype analysis plots for all peptides in the res/biotype. For more information, please refer to the :ref:`normal_mode_example`.

	Translation mode:
	

	In translation mode, BamQuery will expect to find in the input folder `path_to_input_folder` the files BAM_directories.tsv, peptides.tsv and BAM_ribo_directories.tsv. BamQuery will output, in addition to the transcript expression level (RNA bam files), the translation level (ribosome profile files). In this mode, BamQuery can be used as a means to verify the translation of peptides of interest. To do this, BamQuery will search for peptide locations in the BAM/CRAM file(s) in the BAM_directories.tsv and also in the BAM_ribo_directories.tsv directories. Along with the expression heatmaps of each peptide for each BAM/CRAM file/files, you will find the biotype analysis plots for all peptides in the res/biotype. For more information, please refer to the translation mode example.


--strandedness
--------------

When using this option, BamQuery will take into account the strand on which the peptide locations are located. For this, BamQuery will take into account the strandability of each bam file to count reads according to the strand of the queried genomic positions. This takes into account the strandedness of the bam files, so the library (stranded/non-stranded, pair-end, single-end, forward or reverse direction) is automatically detected for each bam file.

If the option is not included, all bam files will be treated according to the pair-end, single-end library but in unstranded mode.


--th_out
--------

The th_out option changes the threshold that is considered for comparing the expression levels of different tissues. By default, this threshold is 8.55 rphm (reads per hundred million). 

--light
-------

In this mode, BamQuery will only display the peptide count and normalization. Therefore, no biotyping analysis will be performed for the peptides. For more information, see the :ref:`light_mode_example`.

--dbSNP
-------

This option allows to choose between three versions of dbSNPs: 149 / 151 / 155. dbSNP 149 is the default. If you don't want to use any release specify 0 for this argument.

--c
---
This option allows to choose between the most COMMON SNPs from the dbSNP release that you choose with the argument above.

--plots
-------
This option sets BamQuery to produce pie charts in the biotype analysis step.

--genome_version
----------------
This option allows to choose between three genome versions : v26_88 / v33_99 / v38_104. genome version v26_88 is the default. 




