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
		  --c                   Take into account the only common SNPs from the dbSNP
		                        database chosen
		  --strandedness        Take into account strandedness of the samples
		  --light               Display only the count and norm count for peptides and
		                        regions
		  --sc                  Query Single Cell Bam Files
		  --var                 Keep Variants Alignments
		  --maxmm               Allows STAR to generate high amount of alignments
		  --overlap             Count overlapping reads
		  --plots               Plot biotype pie-charts
		  --dev                 Save all temps files

====================



Positional Arguments
====================


**path_to_input_folder**
-------------------------

To run BamQuery, create an input folder that includes two files: BAM_directories.tsv and peptides.tsv. 

.. code::

	
	Input
	----|
	    |---:> BAM_directories.tsv
	    |---:> peptides.tsv
	    

The BAM_directories.tsv file must contain the BAM/CRAM files in which the peptides are queried. |br|

The peptides.tsv file must contain the peptides for which BamQuery performs the RNA expression search in the BAM/CRAM files.
(See: :ref:`format input files`).


**name_exp**
-------------

Name of the search as a means of identification.

----------------


Optional Arguments
==================


**-\-mode**
------------

BamQuery has to modes of search : normal or translation. |br|
By default, BamQuery runs in normal mode.

	**Normal mode:**
	In normal mode, BamQuery expects to find in BAM_directories.tsv the BAM/CRAM files corresponding to the RNA-seq datasets. 
	For more information, see :ref:`normal mode example`.

	**Translation mode:**
	Instead of BAM_directories.tsv, BamQuery expects to find a BAM_ribo_directories.tsv corresponding to the Ribo-seq datasets. In this mode, BamQuery can be used as a means to verify the presence of ribosomal profile reads for the peptides of interest. 
	For more information, see :ref:`translation_mode_example`.


.. _genome version:

**-\-genome_version**
----------------------
This option allows to choose between three genome versions supported by BamQuery : v26_88 / v33_99 / v38_104. |br|
By default, genome version v26_88. 


**-\-th_out**
--------------

The th_out option changes the threshold that is considered to highlight in black the heatmap boxes representing RNA expression in Bam files where a peptide has an average rphm>th_out. |br|
By default, this threshold is 8.55 rphm (reads per hundred million). 

.. _dbsnp:

**-\-dbSNP**
-------------

This option allows you to choose between three versions of dbSNPs: 149 / 151 / 155. To specify that no dbSNP version shoubl be used, use dbSNP=0. |br|
By default, dbSNP 149. 


**-\-c**
---------
This option allows only to choose the most COMMON SNPs from the dbSNP release that you choose with the argument above.


.. _strandedness:

**-\-strandedness**
--------------------

When using this option, BamQuery takes into account the strand on which the peptide is located in the genomic location to count the overlapping reads. 
For each Bam file, BamQuery automatically detects the library (stranded/non-stranded, pair-end, single-end, forward or reverse direction). |br|
By defatul, all bam files will be treated according to the pair-end, single-end library but in unstranded mode.


**-\-light**
-------------

In this mode, BamQuery only displays peptide counting and normalization. Therefore, no biotyping analysis will be performed for peptides. |br| 
For more information, see :ref:`light_mode_example`.

**-\-sc**
---------

BamQuery expects to find in BAM_directories.tsv the BAM/CRAM files corresponding to the single cell RNA-seq datasets. BamQuery reports the expression of each peptide in cell populations and generates specific output. |br| 
For more information, see :ref:`single_cell_example`.


**-\-var**
----------
This option sets BamQuery to keep variant alignments where the genome reference translates exactly for the peptide even if the aligned MCS contains mismatches and are not supported by any annotated SNPs. |br| 
For more information, see :ref:`variant_aligments`.

**-\-maxmm**
------------
This option changes some of the STAR parameters (in the MCS alignment process, see :ref:`collect locations`) to allow STAR to generate a large number of alignments. |br|
The new values for the modified STAR parameters are: |br|

.. code::

	--winAnchorMultimapNmax 100000
	--outFilterMultimapNmax 100000
	--limitOutSAMoneReadBytes 26600000
	--outFilterMultimapScoreRange 2
	--alignTranscriptsPerReadNmax 100000


.. warning::
	With this option the STAR aligner will take longer to align the MCS with the genome.


**-\-overlap**
--------------
BamQuery counts an RNA-seq read if the read completely spans the MCS, however, with this option BamQuery also counts RNA-seq reads that overlap at least 60% of the MCS. 


**-\-plots**
-------------
This option sets BamQuery to produce pie charts in the biotype analysis step.


**-\-dev**
----------
This option allows you to save all intermediate files.

.. warning::
	Intermediate files can take up a lot of space.



.. |br| raw:: html

      <br>


