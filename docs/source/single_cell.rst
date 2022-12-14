********************
single cell example
********************

.. _single_cell_example:
BamQuery single cell was designed to perform searches in BAM files from single-cell RNA-seq data. Therefore, BamQuery only reports peptide counts for each cell in each single-cell BAM file, so no biotyping analysis is performed for peptides (no graph (heat map, pie chart) is produced).


**Input**
#########



**Command line:**

.. code::

	BamQuery.py ./sc_example/Input sc_example --sc

As for normal mode, the input folder `path_to_input_folder` must containt the files : **BAM_directories.tsv** and **peptides.tsv**.

**Output**
##########

BamQuery creates an **output** directory in the same path as the input folder.

This directory contains 3 folders and the main results are organized as follows:

.. code::

	├── alignments
	│   ├── missed_peptides.info
	│   └── sc_example_info_alignments.xlsx
	├── logs
	│  ├── BamQuery_Res_sc_example.log
	│  └── Information_BAM_directories.log
	└── res
	    ├── sc_example_rna_sc_count_All_alignments.csv
	    └── sc_example_rna_sc_count.csv

The output files with `BamQuery single cell` are similar to those in `normal mode`, see `output_normal_mode_example`_ for detailed information about the output files : missed_peptides.info, sc_example_info_alignments.xlsx and logs.

-------------

**res**
=======

.. code::

	res
	    ├── sc_example_rna_sc_count_All_alignments.csv
	    └── sc_example_rna_sc_count.csv

.. _sc_example_rna_sc_count_All_alignments:

`sc_example_rna_sc_count_All_alignments.csv`: 
Reports for each MCS of each peptide at a given location the total of RNA-seq reads in each cell found in single-cell BAM file included in **BAM_directories.tsv**. |br| 
Each BAM file is identified in sc_example_rna_sc_count_All_alignments.csv according to the order in which the BAM file was included in the **BAM_directories.tsv**. For instance: `0_TCTGAGACAGGTCGTC` means that cell `TCTGAGACAGGTCGTC` is included in the first (`0`) BAM file in **BAM_directories.tsv**.

.. thumbnail:: _images/sc_example_rna_sc_count_All_alignments.jpg

.. _sc_example_rna_sc_count:

`sc_example_rna_sc_count.csv`: 
Reports for each peptide the total of RNA-seq reads in each cell found in single-cell BAM file included in **BAM_directories.tsv**. |br| 
Each BAM file is identified in sc_example_rna_sc_count_All_alignments.csv according to the order in which the BAM file was included in the **BAM_directories.tsv**. For instance: `0_TCTGAGACAGGTCGTC` means that cell `TCTGAGACAGGTCGTC` is included in the first (`0`) BAM file in **BAM_directories.tsv**.

.. thumbnail:: _images/sc_example_rna_sc_count.jpg

