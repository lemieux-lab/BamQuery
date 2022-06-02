========
Methods
========

.. _bamquery_steps:

BamQuery Steps
==============

BamQuery follows four important steps for each peptide queried. 

.. image:: _static/Approach.png
   :alt: Approach
   :align: center


1. Collect the locations of the MAP Coding Sequences (MCS) in the genome. 
------------------------------------------------------------------------

The sequences from the manual mode peptide reverse translation will be compiled into a **<query_name>.fastq** file along with any MCS in MAP coding sequence mode. These sequences will then be mapped to the reference genome version (GENCODE version 26, 33 or 38 ) selected by the user using the STAR aligner.
**MCS genomic locations** are collected from the perfect alignments in the output STAR file Aligned.out.sam. 
A perfect alignment is then defined as one in which the MCS exactly matches the sequence at the alignment according to the reference genome. A perfect alignment is also one in which the MCS does not exactly match the reference, but any mismatches are supported by SNPs included in the user-selected dbSNP (149,151,155 dbSNP releases). 


2. Retrieve the total count of RNA-seq reads for the MAP. 
---------------------------------------------------------

Each BAM/CRAM file will be queried at all genomic locations collected for a MAP. 
Primary RNA-seq reads at each location will be queried using the python **pysam** library. Each read will then be scrutinized to count only those reads that carry the MCS reported at the location.
Finally, BamQuery will report the total count of the primary RNA-seq reads on each BAM/CRAM file that exactly match the MCSs of a peptide. 


.. warning::
	Total RNA-seq reads count changes according to the strandedness parameter.

3. Normalization. 
-----------------

The :math:`tr_{MAP}` count is transformed into a number of reads detected per :math:`10^{8}` reads sequenced (𝑟𝑝ℎ𝑚) following the formula : :eq:`rphm`. :math:`R_{t}` represents the total number of reads sequenced in a given RNA-Seq dataset. These final values are log-transformed :math:`log_{10} (𝑟𝑝ℎ𝑚 + 1)`.


.. math::
	𝑟𝑝ℎ𝑚 = \frac{tr_{MAP} } {R_{t} } * 10^{8} 
	:label: rphm

4. Biotype classification. 
--------------------------

Each MAP is classified according to all its MCS genomic locations and their transcription level. Biotype screening is performed using GENCODE annotations according to the reference genome selected by the user and the RepeatMasker annotations to account for human ERE sequences. All MCS genomic locations are included in a BED file for intersect with annotations using BEDtools 73.  To account for locations overlapping several transcripts, reads were scales based on 


Format Input Files
===================

	Let's take a look at the format of these files.

**a. BAM_directories.tsv**
--------------------------

	This file should look like follows:

	.. image:: _images/BAM_directories.png
	   :alt: Format BAM_directories.tsv
	   :align: left

	BamQuery will search for all BAM/CRAM files in each path you include.

	Note that:

	1. The first column is the name of the BAM/CRAM file or group of files to be queried. This name should describe the type of BAM/CRAM file(s).
	2. The second column should be the path to the BAM/CRAM file(s).
	3. The first and second columns are separated by a tab space. Please make sure that the space between the two columns is a tab space, otherwise BamQuery will throw an exception.
	4. 4. Do not use any headers in your tsv file.


**b. peptides.tsv**
-------------------

	This file should look like follows:

	.. image:: _images/peptides_file_format.png
	   :alt: Format peptides.tsv
	   :align: left


	Note that all peptides can be pulled into a single peptides.tsv, however, you must follow the format assigned for each mode.

	**Peptides in peptide mode:** 
	running peptides in this mode assumes that you only know the amino acid sequence of the peptide. BamQuery will reverse translate those peptides in order to find all possible locations in the genome for the MAP coding sequences.
	
	Two columns separated by a tab space: 

	a. In the first column add the amino acid sequence of the peptide.
	b. In the second column add the type of peptide to identify it. This name, for example, may refer to the condition or sample in which the peptide was identified. 
		
	**Peptides in coding sequence mode:** running peptides in this mode assumes that the amino acid and nucleotide sequence of each peptide queried is known.
	
	Three columns separated by a tab space: 

	a. In the first column the amino acid sequence of the peptide is added. 
	b. In the second column the nucleotide sequence of the peptide.
	c. In the third column the type of peptide to identify it. 
		
	**Peptides in manual mode:** 
	running peptides in this mode assumes that the amino acid, nucleotide coding sequence, and the location and strand of each peptide queried are known.
	
	Five columns separated by a tab space:

	a. In the first column add the amino acid sequence of the peptide.
	b. In the second column add the nucleotide sequence of the peptide.
	c. In the third column add the position of the peptide.
	d. In the fourth column add the strand (-) backward or (+) forward for the location of the peptide in the genome.
	e. In the last column add the type of peptide to identify it. 

	.. warning::
		The peptide location must follow the format: chrX:start-end|start-end. Note: chrX (for any chromosome), start = start location, end = end location. Only use "|" to specify if the peptide is spliced.
		The strand must be specified as (-) backward or (+) forward.
		


.. warning::
	If a peptide has several peptide types, separate each peptide type with ", or ;". For example: `lymphoma,colon`, would mean that the peptide was identified in lymphoma and colon cells.

