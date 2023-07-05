
********************
single cell example
********************

-------

BamQuery single cell was designed to perform searches in BAM files from single-cell RNA-seq data. Therefore, BamQuery only reports peptide counts for each cell in each single-cell BAM file, so no normalization (rphm values) and biotyping analysis is performed for peptides (no graph (heat map, pie chart) is produced).


**Input**
#########


**Command line:**

.. code::

	BamQuery.py ./sc_example/Input sc_example v38_104 --sc

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
	    ├── sc_example_cell_names_identification.csv
		├── sc_example_index_bam_files.csv
		├── sc_example_rna_sc_count_All_alignments.csv
	    └── sc_example_rna_sc_count.csv

The output files with `BamQuery single cell` are similar to those in `normal mode`, see `output_normal_mode_example`_ for detailed information about the output files : missed_peptides.info, sc_example_info_alignments.xlsx and logs.

-------------

:purple:`res`
=============

.. code::

	res
	    ├── sc_example_cell_names_identification.csv
	    ├── sc_example_index_bam_files.csv
	    ├── sc_example_rna_sc_count_All_alignments.csv
	    └── sc_example_rna_sc_count.csv

.. _sc_example_cell_names_identification:

**sc_example_cell_names_identification.csv**: |br| 
Reports the cells in each single-cell BAM file included in **BAM_directories.tsv**. |br| 
Each BAM file is associated with an index, and each cell is associated with an ``id``. 

.. thumbnail:: _images/sc_example_cell_names_identification.jpg

.. _sc_example_index_bam_files:

**sc_example_index_bam_files.csv**: |br| 
Reports the index associated to each single-cell BAM file included in **BAM_directories.tsv**. |br| 

.. thumbnail:: _images/sc_example_index_bam_files.jpg

.. _sc_example_rna_sc_count_All_alignments:

**sc_example_rna_sc_count_All_alignments.csv**: |br| 
Reports for each MCS at a given location the total of RNA-seq reads in each cell found in single-cell BAM file included in **BAM_directories.tsv**. |br| 
Each cell is identified by its assigned ``id``. |br| 
For instance: ``1`` corresponds to the cell ``AAACCTGAGACACGAC`` which is included in the first (``0``) BAM file in **BAM_directories.tsv** ``HCATisStab7509734_outs``.

.. thumbnail:: _images/sc_example_rna_sc_count_All_alignments.jpg

.. _sc_example_rna_sc_count:

**sc_example_rna_sc_count.csv**: |br| 
Reports for each peptide the total of RNA-seq reads in each cell found in single-cell BAM file included in **BAM_directories.tsv**. |br| 
Each cell is identified by its assigned ``id``. |br| 
For instance: ``1`` corresponds to the cell ``AAACCTGAGACACGAC`` which is included in the first (``0``) BAM file in **BAM_directories.tsv** ``HCATisStab7509734_outs``.

.. thumbnail:: _images/sc_example_rna_sc_count.jpg

.. note::
	Use BamQuery to collect the unique molecular identifier (UMI) count for each peptide by using the following command:

	.. code::

		BamQuery.py ./sc_example/Input sc_example v38_104 --sc --umi


