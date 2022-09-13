========
Methods
========

.. _bamquery_steps:

BamQuery Steps
==============

BamQuery works on Bam files in five steps. 

.. image:: _images/Approach.png
   :alt: Approach
   :align: center


1. Reverse translation of MHC-I-associated peptides (MAPs). 
-----------------------------------------------------------

Each MAP in `Peptide mode`_ peptide mode is reverse-translated. MAP coding sequences (MCS) are compiled into a **<query_name>.fastq** file along with any MCS from MAPs in the `MAP coding sequence mode`_ .

.. _collect locations:

2. Collect the locations of the MAP Coding Sequences (MCS) in the genome. 
--------------------------------------------------------------------------

The MCSs in the **<query_name>.fastq** file are mapped through the STAR aligner to the reference genome version (GENCODE version 26, 33 or 38 ) selected by the user.

**The genomic locations of the MCS** are collected from the perfect alignments in the output STAR file Aligned.out.sam. 
A perfect alignment is defined as one in which the MCS exactly matches the reference genome sequence. A perfect alignment is also one in which the MCS does not exactly match the reference, but any mismatch is supported by SNPs included in the user-selected dbSNP (149,151,155 dbSNP releases). 


3. Count RNA-seq reads that exactly bear the MCS. 
--------------------------------------------------

Each BAM/CRAM file will be queried at all genomic locations collected for a MAP. 
The primary RNA-seq reads at each location will be queried using the python **pysam** library. Each read will be examined to count only those reads that carry the MCS reported at the location.
Finally, BamQuery will report the total count of RNA-seq primary reads in each BAM/CRAM file that exactly match the MCS of a peptide. 


.. warning::
	Total RNA-seq reads count changes according to the :ref:`strandedness` parameter.

4. Normalization. 
-----------------

The :math:`tr_{MAP}` count is transformed into a number of reads detected per :math:`10^{8}` reads sequenced (𝑟𝑝ℎ𝑚) following the formula : :eq:`rphm`. :math:`R_{t}` represents the total number of reads sequenced in a given RNA-Seq dataset. These final values are log-transformed :math:`log_{10} (𝑟𝑝ℎ𝑚 + 1)` to allow comparison and averaging between samples, thus removing the bias of large values.


.. math::
	𝑟𝑝ℎ𝑚 = \frac{tr_{MAP} } {R_{t} } * 10^{8} 
	:label: rphm


.. _biotype:

5. Biotype classification. 
--------------------------

Each MAP is classified according to all its MCS genomic locations and their transcription level. 
Biotype screening is performed using GENCODE annotations according to the reference genome selected by the user and the RepeatMasker annotations to account for human ERE sequences. 
All MCS genomic locations are included in a BED file for intersect with annotations using BEDtools.
To account for locations overlapping with several transcripts, reads were scaled according to the read distribution coefficient for each biotype from estimation of those parameters with the expectation maximization (EM) statistical model. |br|
For more information, see :ref:`biotypes`.


As a result, BamQuery attributes a comprehensive RNA expression to any MAP of interest in any user-selected RNA-seq dataset. 


---------------


Input Modes
===================

BamQuery is designed to analyze MAPs ranging in length from 8 to 11 amino acids (aa). 
As peptide input, BamQuery supports three different formats that can be pulled into a single input file (See `peptides tsv`_). 

.. _Peptide mode:

**A) Peptide mode:** only the amino acid sequence of the MAP is provided, hence BamQuery performs a comprehensive search for its RNA-seq expression. 

.. _MAP coding sequence mode:

**B) MAP coding sequence mode:** the amino acid sequence of the MAP is provided, hence BamQuery performs the search for the expression of the given MCS. 

**C) Manual mode:** the amino acid sequence of the MAP is provided followed by a MCS, the corresponding location in the genome of the given MCS, the strand (+ forward or - reverse), whereby BamQuery performs the expression search at the given location for the given MCS.



-----------


.. _format input files:

Format Input Files
===================


BamQuery requires two input file paths to search for RNA expression:

**A) BAM_directories:** list of Bam files in which the search is performed
**B) peptides.tsv:** list of peptides to be searched


**A) BAM_directories.tsv**
--------------------------

	This file should look like follows:

	.. image:: _images/BAM_directories.png
	   :alt: Format BAM_directories.tsv
	   :align: left

	BamQuery collects all BAM/CRAM files in each path included in the list. For instance from the path /home/gtex/, BamQuery collects all the bam files for every tissue in gtex.

	Note that:

	1. The first column is the name of the BAM/CRAM file or group of files to be queried. This name should describe the type of BAM/CRAM file(s).
	2. The second column should be the path to the BAM/CRAM file(s).
	3. The first and second columns are separated by a tab space. 
	4. Do not use any headers in your tsv file.


.. _peptides tsv:


**B) peptides.tsv**
-------------------

	This file should look like follows:

	.. image:: _images/peptides_file_format.png
	   :alt: Format peptides.tsv
	   :align: left


	Note that all modes can be merged into a single peptides.tsv, however, you must follow the format assigned for each mode.

	
	.. warning::
		If a peptide has several peptide types, separate each peptide type with ", or ;". For example: `lymphoma,colon`, would mean that the peptide was identified in lymphoma and colon cells.

	**Peptides in peptide mode:** |br|
	Two columns separated by a tab space: |br|
	a. amino acid sequence of the peptide. |br|
	b. type of peptide to identify it. This name, for example, may refer to the condition or sample in which the peptide was identified. 
		
	**Peptides in coding sequence mode:** |br|
	Three columns separated by a tab space: |br|
	a. amino acid sequence of the peptide. |br|
	b. nucleotide sequence of the peptide. |br|
	c. type of peptide to identify it. 
		
	**Peptides in manual mode:** |br|
	Five columns separated by a tab space: |br|
	a. amino acid sequence of the peptide. |br|
	b. nucleotide sequence of the peptide. |br|
	c. position of the peptide. |br|
	d. strand backward (-) or forward (+) for the location of the peptide in the genome. |br|
	e. type of peptide to identify it. 

	.. warning::
		The peptide location must follow the format: chrX:start-end|start-end. Note: chrX (for any chromosome), start = start location, end = end location. Only use "|" to specify if the peptide is spliced.
		The strand must be specified as (-) backward or (+) forward.
		

.. |br| raw:: html

      <br>