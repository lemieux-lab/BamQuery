====================
Command Line Options
====================

At the command line::

    BamQuery --help


.. code::

		usage: BamQuery.py [-h] [--mode MODE] [--th_out TH_OUT] [--dbSNP DBSNP] [--c]
						[--strandedness] [--light] [--sc] [--var] [--maxmm]
						[--overlap] [--plots] [--m] [--dev] [--t T]
						path_to_input_folder name_exp genome_version

		======== BamQuery ========

		positional arguments:
			path_to_input_folder  Path to the input folder where to find BAM_directories.tsv and peptides.tsv
			name_exp              BamQuery search Id
			genome_version        Genome human releases : v26_88 / v33_99 / v38_104; Genome mouse releases : M24 / M30

		optional arguments:
						-h, --help            show this help message and exit
						--mode MODE           BamQuery search mode : normal / translation
						--th_out TH_OUT       Threshold to assess expression comparation with other tissues
						--dbSNP DBSNP         Human dbSNP : 149 / 151 / 155
						--c                   Take into account the only common SNPs from the dbSNP database chosen
						--strandedness        Take into account strandedness of the samples
						--light               Display only the count and norm count for peptides and regions
						--sc                  Query Single Cell Bam Files
						--umi                 Count UMIs in Single Cell Bam Files
						--var                 Keep Variants Alignments
						--maxmm               Enable STAR to generate a larger number of alignments
						--overlap             Count overlapping reads
						--plots               Plot biotype pie-charts
						--m                   Mouse genome
						--dev                 Save all temps files
						--t T                 Specify the number of processing threads to run BamQuery. The default is 4


====================



Positional Arguments
====================


**path_to_input_folder**
-------------------------

To run BamQuery, create an input folder that includes two files: BAM_directories.tsv and peptides.tsv. 

.. code::

	
	├── Input
		├── BAM_directories.tsv
		└── peptides.tsv
	    

The BAM_directories.tsv file must contain the BAM/CRAM files in which the peptides are queried. |br|

The peptides.tsv file must contain the peptides for which BamQuery performs the RNA expression search in the BAM/CRAM files.
(See: :ref:`format input files`).


**name_exp**
-------------

Name of the search as a means of identification.

.. _genome version:

**genome_version**
-------------------

This option allows to choose between three GRCh38 human genome annotation publication versions: v26_88 / v33_99 / v38_104. |br|

This option together with **-\-m** (which allows BamQuery to search in the mouse genome) supports two mouse genome annotation versions of GRCm38 and GRCm39, respectively: M24 / M30. |br|

You need to donwload any of the human or mouse supported versions in BamQuery. |br|
Genome human releases : v26_88 / v33_99 / v38_104 |br|
Genome mouse releases : M24 / M30 |br|

Please see :ref:`installation` to follow the installation instructions for the genome versions..

.. note::
	**Gencode Human Releases supported with BamQuery**

	Gencode_human/release_26 : v26_88 |br|
	Evidence-based annotation of the human genome (GRCh38), version 26 (Ensembl 88)
	https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/


	Gencode_human/release_33 : v33_99 |br|
	Evidence-based annotation of the human genome (GRCh38), version 33 (Ensembl 99)
	https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/


	Gencode_human/release_38 : v38_104 |br|
	Evidence-based annotation of the human genome (GRCh38), version 38 (Ensembl 104)
	https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/



	**Gencode Mouse Releases supported with BamQuery**

	Gencode_mouse/release_M24 : M24 |br|
	Evidence-based annotation of the mouse genome (GRCm38), version M24 (Ensembl 99)
	https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M24/


	Gencode_mouse/release_M30 : M30 |br|
	Evidence-based annotation of the mouse genome (GRCm39), version M30 (Ensembl 107)
	https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M30/

----------------


Optional Arguments
==================



**-\-mode**
------------

BamQuery has two modes of search : normal or translation. |br|
By default, BamQuery runs in normal mode.

	**Normal mode:**
	In normal mode, BamQuery expects to find in BAM_directories.tsv the BAM/CRAM files corresponding to the RNA-seq datasets. 
	For more information, see :ref:`normal_mode_example`.

	**Translation mode:**
	Instead of BAM_directories.tsv BamQuery expects a BAM_ribo_directories.tsv corresponding to the Ribo-seq datasets. In this mode, BamQuery can be used as a means to verify the presence of ribosomal profiling reads for the peptides of interest. 
	For more information, see :ref:`translation_mode_example`.


**-\-th_out**
--------------

The th_out option changes the threshold that is considered to highlight in black the heatmap boxes representing RNA expression in Bam files where a peptide has an average rphm>th_out. |br|
By default, this threshold is 8.55 rphm (reads per hundred million). 

.. _dbsnp:

**-\-dbSNP**
-------------

This option allows you to choose between three versions of dbSNPs: 149 / 151 / 155. |br|
By default, dbSNP 0. 


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

BamQuery expects to find in BAM_directories.tsv the BAM/CRAM files corresponding to the single cell RNA-seq datasets. |br|
BamQuery reports the RNA-seq read count for each peptide in cell populations and generates specific output. |br| 
For more information, see :ref:`single_cell_example`.

**-\-umi**
-----------

BamQuery expects to find in BAM_directories.tsv the BAM/CRAM files corresponding to the single cell RNA-seq datasets. |br|
BamQuery reports instead of RNA-seq read count, it reports the unique molecular identifier (UMI) count for each peptide in cell populations and generates specific output. |br| 
For more information, see :ref:`single_cell_example`.

**-\-var**
----------
A variant alignment refers to an alignment where the mapped MCS may deviate from the reference genome sequence by a maximum of 3 nucleotides (1 amino acid). |br| 
In these cases, single nucleotide variants are taken into account even though they are not included in the selected dbSNP.  

.. note::
	Allowing BamQuery to maintain variant alignments could facilitate the evaluation of the expression of mutated MAPs. 
	However, using this option generates a large number of alignments that would impact execution time.

**-\-maxmm**
------------
This option changes some of the STAR parameters (in the MCS alignment process, see :ref:`collect locations`) to allow STAR to generate a larger number of alignments. |br|
The new values for the modified STAR parameters are: |br|

.. code::

	--winAnchorMultimapNmax 20000
	--outFilterMultimapNmax 100000
	--outFilterMultimapScoreRange 4
	--alignTranscriptsPerReadNmax 100000
	--seedPerWindowNmax 1500
	--seedNoneLociPerWindow 1500
	--alignWindowsPerReadNmax 20000
	--alignTranscriptsPerWindowNmax 1500

.. warning::
	With this option the STAR aligner will take longer to align the MCS with the genome.


**-\-overlap**
--------------
BamQuery counts an RNA-seq read if the read completely spans the MCS, however, with this option BamQuery also counts RNA-seq reads that overlap at least 60% of the MCS. 


**-\-plots**
-------------
This option sets BamQuery to produce pie charts in the biotype analysis step.

**-\-m**
-------------
This option sets BamQuery to search for peptides in the mouse genome. |br|
Along with the **-\-genome_version** option BamQuery can be parameterized to run the search on either of the two supported GRCm38 and GRCm39 mouse genome annotation versions: M24 / M30. If only the **-\-m** option is passed as an argument, BamQuery takes the default M24 mouse genome annotation version. |br|
By default, the mouse genome annotation versions: M24 / M30, are used with the EVA database of genomic variation for the GRCm38 and GRCm39, respectively.


**-\-dev**
----------
This option allows you to save all intermediate files.

.. warning::
	Intermediate files can take up a lot of space.


**-\-t**
----------
Specify the number of processing threads to run BamQuery. By default, BamQuery runs in 4 threads.


.. |br| raw:: html

      <br>


