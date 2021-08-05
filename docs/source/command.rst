====================
Command Line Options
====================

At the command line::

    BamQuery.py --help


.. code-block:: bash

	    usage: BamQuery.py [-h] [--mode MODE] [--strandedness] [--th_out TH_OUT]
	                   [--light]
	                   path_to_input_folder name_exp

		======== BamQuery ========

		positional arguments:
		  path_to_input_folder  Path to the input folder where to find
		                        BAM_directories.tsv and peptides.tsv
		  name_exp              BamQuery search Id

		optional arguments:
		  -h, --help            show this help message and exit
		  --mode MODE           BamQuery search mode : normal / translation
		  --strandedness        Take into account strandedness of the samples
		  --th_out TH_OUT       Threshold to assess expression comparation with other
		                        tissues
		  --light               Display only the count and norm count for peptides and
		                        regions

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

====================

Format Input Files
===================

	Let's take a look at the format of these files.

a. BAM_directories.tsv
----------------------

	This file should look like follows:

	.. image:: _images/Bam_directories.png
	   :alt: Format BAM_directories.tsv
	   :align: left

	BamQuery will search for all BAM/CRAM files in each path you include.

	Note that:

	1. The first column is the name of the BAM/CRAM files you want to query. The name describes the type of BAM/CRAM file(s)/file(s) found in the path.
	2. The second column is the path to the BAM/CRAM file(s)/file(s).
	3. The first and second columns are separated by a tab space. Please make sure that the space between the two columns is a tab space, otherwise BamQuery will throw an exception.
	4. Do not use any headers in your tsv file.


b. peptides.tsv
---------------

	This file should look like follows:

	.. image:: _images/peptides_tsv.png
	   :alt: Format peptides.tsv
	   :align: left


	Please note that all the peptides can be pull into a single peptides.tsv, howver, you must follow the assigned format for each mode.

	1. Peptides in peptide mode: running peptides in this mode assumes that you only know the amino acid sequence of the peptide. BamQuery will search for you all the possible coding sequences of every peptide and find all possible positions in the genome for the coding sequences.
		a. Two columns separated by a tab space: in the first column add the amino acid sequence of the peptide, in the second column add the type of peptide to identify it. This name, for example, can refer to the condition or sample in which the peptide was identified. 
		b. If a peptide has several peptide types, separate each peptide type with "_". For example: `lymphoma_colon`, would mean that the peptide was identified in lymphoma and colon cells. 
		c. BamQuery will reverse translate all peptides in this mode.

	2. Peptides in coding sequence mode: running peptides in this mode assumes that the amino acid and nucleotide sequence of each peptide being queried is known.
		a. Three columns separated by a tab space: in the first column add the amino acid sequence of the peptide, in the second column add the nucleotide sequence of the peptide, and in the third column add the type of peptide to identify it. This name, for example, can refer to the condition or sample in which the peptide was identified.
		b. If a peptide has several peptide types, separate each peptide type with "_". For example: `lymphoma_colon`, would mean that the peptide was identified in lymphoma and colon cells.

	3. Peptides in manual mode: running peptides in this mode assumes that the amino acid, nucleotide coding sequence, and the location and chain of each peptide being queried are known.
		a. Five columns separated by a tab space: in the first column add the amino acid sequence of the peptide, in the second column add the nucleotide sequence of the peptide, in the third column add the position of the peptide, in the fourth column add the string (-) backward or (+) forward for the location of the peptide in the genome, and in the last column add the type of peptide to identify it. This name, for example, can refer to the condition or sample in which the peptide was identified.
		b. The peptide position has the following format: chrX:start-end|start-end. Note: chrX (for any chromosome), start = start position, end = end position. Only use "|" to specify if the peptide is spliced.
		c. The strand should be specified as follows (-) backward or (+) forward.
		d. If a peptide has multiple peptide types, separate each peptide type with "_". For example: `lymphoma_colon`, would mean that the peptide was identified in lymphoma and colon cells.


Optional Arguments
==================

--mode
------

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

