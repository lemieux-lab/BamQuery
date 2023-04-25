########
Cookbook
########

See this guide for examples of different uses of BamQuery.

BamQuery can search for peptide expression in the human (GRCh38, annotation versions v26_88 / v33_99 / v38_104) 
and mouse (GRCm38 and GRCm39, annotation versions M24, M30, respectively) genomes.


=======================


.. _normal_mode_example:

*******************
normal mode example
*******************

**Input**
#########

**Command line:**

.. code::

	BamQuery.py ./normal_mode_example/Input normal_mode_example v38_104

Input folder `path_to_input_folder` containing the files : **BAM_directories.tsv** and **peptides.tsv**. 
Please download BAM_directories.tsv and peptides.tsv to access an example of the format of these files.

See:
:download:`BAM_directories.tsv <_static/BAM_directories.tsv>`

See:
:download:`peptides.tsv <_static/peptides.tsv>`

.. _output_normal_mode_example:

**Output**
##########

BamQuery creates an **output** directory in the same path as the input folder.

This directory contains 4 folders and the main results are organized as follows:

.. code::

	├── alignments
	│   ├── missed_peptides.info
	│   └── normal_mode_info_alignments.xlsx
	├── logs
	│   ├── BamQuery_Res_normal_mode.log
	│   └── Information_BAM_directories.log
	├── plots
	│   ├── biotypes
	│   │   ├── biotype_by_sample_group
	│   │   │   ├── all_peptides
	│   │   │   │   └── normal_mode_All_peptides.pdf
	│   │   │   └── by_peptide_type
	│   │   │       ├── normal_mode_0D5P_Melanoma_ctrl_All_samples.pdf
	│   │   │       └── normal_mode_0D5P_Melanoma_DAC_All_samples.pdf
	│   │   └── genome_and_ERE_annotation
	│   │       ├── all_peptides
	│   │       │   └── normal_mode_All_peptides.pdf
	│   │       └── by_peptide_type
	│   │           ├── normal_mode_0D5P_Melanoma_ctrl.pdf
	│   │           └── normal_mode_0D5P_Melanoma_DAC.pdf
	│   └── heat_maps
	│       └── transcription_evidence_heatmap
	│           ├── average_transcription_expression_heatmap
	│           │   ├── normal_mode_rna_norm_all_tissues.pdf
	│           │   ├── normal_mode_rna_norm_selected_tissues.pdf
	│           │   └── norm_info.csv
	│           └── total_transcription_expression_heatmap
	└── res
	    ├── biotype_classification
	    │   ├── full_info_biotypes
	    │   │   ├── 1_Genomic_and_ERE_Annotations_Full.csv
	    │   │   ├── 2_Genomic_and_ERE_Annotations_Summary_Full.csv
	    │   │   └── 3_Genomic_and_ERE_Anno_by_Region_Full.csv
	    │   └── summary_info_biotypes
	    │       ├── 1_General_Gen_and_ERE_Biotype_Consensus.csv
	    │       ├── 2_Sample_Gen_and_ERE_Biotype_Consensus.csv
	    │       └── 3_Group_Samples_Gen_and_ERE_Biotype_Consensus.csv
	    └── normal_mode_count_norm_info.xlsx




.. _alignments:


**alignments**
==============

.. code::

	├── alignments
	   ├── missed_peptides.info
	   └── normal_mode_info_alignments.xlsx


**missed_peptides.info file**


This file reports the peptides for which BamQuery has not been able to find locations to evaluate their expression. 

.. _normal mode example info alignments explanation xlsx file:

**normal_mode_example_info_alignments.xlsx file**


.. _COSMIC: https://cancer.sanger.ac.uk/cosmic

This file reports, for each peptide queried, all locations in the genome that are perfect alignments for one or more MAP coding sequences. In addition, it reports the somatic mutations annotated in the `COSMIC`_ database, encountered in the queried peptides.


`Sheet : Perfect Alignments`

For each MAP: position, strand, MCS, reference amino acid, nucleotide differences and SNVs annotated in the dbSNP database are reported. In the case of the latter, it will be reported that the coding sequence has some difference with the reference genome, but has been compensated by one or more annotated SNVs.


.. thumbnail:: _images/normal_mode_example_info_alignments.jpg


`Sheet : COSMIC Information`

Reports whether the SNVs displayed by the peptides have already been annotated as somatic mutations in the COSMIC database.


.. thumbnail:: _images/Cosmic.jpg

Find here the `Additional information`_ description from COSMIC.


.. _Additional information: https://cancer.sanger.ac.uk/cosmic/download


-----

.. _Logs:

**logs**
========

.. code::

	├── logs
	│   ├── BamQuery_Res_normal_mode.log
	│   └── Information_BAM_directories.log
	      		       

**BamQuery_Res_normal_mode_example.log file**

This file reports all steps that have been performed in the BamQuery search. Refer to this file for the query time of all peptide alignments in the bams, the number of perfect peptide alignments, and the summary of the parameters used in the search.

**Get_Read_Count_BAM_directories.log file**

This file reports for each BAM/CRAM file in the **BAM_directories.tsv** the total number of primary read counts.


-----

**plots**
=========

The Plots folder contains the heat map and biotype analysis expression plots for all peptides.
If the --plots parameter is specified, pie charts of the biotype classification are produced. 

.. code::

	├── plots
	│   ├── biotypes
	│   │   ├── biotype_by_sample_group
	│   │   │   ├── all_peptides
	│   │   │   │   └── normal_mode_All_peptides.pdf
	│   │   │   └── by_peptide_type
	│   │   │       ├── normal_mode_0D5P_Melanoma_ctrl_All_samples.pdf
	│   │   │       └── normal_mode_0D5P_Melanoma_DAC_All_samples.pdf
	│   │   └── genome_and_ERE_annotation
	│   │       ├── all_peptides
	│   │       │   └── normal_mode_All_peptides.pdf
	│   │       └── by_peptide_type
	│   │           ├── normal_mode_0D5P_Melanoma_ctrl.pdf
	│   │           └── normal_mode_0D5P_Melanoma_DAC.pdf
	│   └── heat_maps
	│       └── transcription_evidence_heatmap
	│           ├── average_transcription_expression_heatmap
	│           │   ├── normal_mode_rna_norm_all_tissues.pdf
	│           │   ├── normal_mode_rna_norm_selected_tissues.pdf
	│           │   └── norm_info.csv
	│           └── total_transcription_expression_heatmap
	

**biotypes**
-------------

The `biotype_by_sample_group` folder contains the biotype assignment based on transcription expression, i.e. the biotype is computed based on those locations where there are underlying RNA-seq reads. For more information please refers to :ref:`biotype` and :ref:`biotypes`

This folder contains pie charts organised as follows:

1) `all_peptides`: Pie chart depicting the general assignment of biotypes for all peptides queried based on their transcription in the BAM/CRAM files consulted.  

.. image:: _images/biotype_transcription.jpg
   :alt: Biotype based on transcription
   :align: center
   :scale: 50 %


2) `by_peptide_type`: Pie charts showing the general assignment of biotypes according to each peptide type (specified in the **peptides.tsv** file) based on transcript expression, i.e. biotype is calculated based on the locations where there are underlying RNA-seq reads


The `genome_and_ERE_annotation` folder contains the biotype assignment regardless of transcript expression, i.e. the biotype assignment for each peptide is calculated based on all locations in the genome (expressed or not).

This folder contains pie charts organised as follows:

1) `all_peptides`: Pie chart depicting the general assignment of biotypes for all peptides queried based on all the locations for all the peptides.  

.. image:: _images/biotype_locations.jpg
   :alt: Biotype based on locations
   :align: center
   :scale: 50 %

2) `by_peptide_type`: Pie charts showing the general assignment of biotypes according to each peptide type (specified in the **peptides.tsv** file) based on all the locations for all the peptides. 


.. _heat maps folder:

**heat_maps**
-------------

This folder contains the heat maps representing the transcription expression levels of all peptides queried.

`average_transcription_expression_heatmap` folder: 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Heat maps representing the mean transcription expression for each peptide queried in the tissues associated with the BAM/CRAM files.

`_norm_all_tissues.pdf` : Heat map representing the average level of transcription expression associated with tissue types, computed from all samples in the tissue.

.. thumbnail:: _images/average_transcription_expression_heatmap.jpg


`_norm_selected_tissues.pdf` : Heat map representing the average level of transcription expression associated with tissue types, computed from selected tissues (short list of tissues).

.. thumbnail:: _images/average_transcription_expression_heatmap_selected.jpg


`norm_info.csv`: reports, for each peptide queried, the mean and median values of rphm in the tissues associated with the BAM/CRAM files.

.. thumbnail:: _images/norm_info_.jpg


`total_transcription_expression_heatmap` folder: 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Heat maps representing the mean transcription expression and total number of RNA-seq reads for each peptide queried in each BAM/CRAM file.

.. warning::
	1. Heat maps are produced for searches with less than 400 peptides.
	2. Heat maps in total_transcription_expression_heatmap is produced only if the number of BAM/CRAM files queried are less than 100.


-----------

**res**
=======

.. code::

	res
	    ├── biotype_classification
	    │   ├── full_info_biotypes
	    │   │   ├── 1_Genomic_and_ERE_Annotations_Full.csv
	    │   │   ├── 2_Genomic_and_ERE_Annotations_Summary_Full.csv
	    │   │   └── 3_Genomic_and_ERE_Anno_by_Region_Full.csv
	    │   └── summary_info_biotypes
	    │       ├── 1_General_Gen_and_ERE_Biotype_Consensus.csv
	    │       ├── 2_Sample_Gen_and_ERE_Biotype_Consensus.csv
	    │       └── 3_Group_Samples_Gen_and_ERE_Biotype_Consensus.csv
	    └── normal_mode_count_norm_info.xlsx


.. _biotype_classification:

**biotype_classification**
--------------------------

.. _Ensembl: https://m.ensembl.org/info/genome/genebuild/biotypes.html

.. note::
	The biotype annotation is derived from the intersection of the peptide genomic locations with Ensembl and ERE annotations. For more information see :ref:`biotypes`.

	From Ensembl annotations, 3 levels of biotypes are reported : gene level, transcript level and genomic position level. 

	At the genetic level, the biotype assigned to the genomic location is given by the biotype of the gene in `Ensembl`_ that overlaps with the location, for instance:
		* protein_coding,
		* lincRNA,
		* intergenic...

	At the transcript level, the biotype assigned to the genomic location is given by the biotype of the transcript in `Ensembl`_ that overlaps with the location, for instance:
		* protein_coding,
		* processed_transcript, TEC, etc...

	At the genomic position level, the biotype assigned to the genomic location is given by the exact location where the peptide is located in the transcript in `Ensembl`_, for example:
		* In_frame,
		* junctions,
		* introns,
		* 3'UTR, etc...

	As for the ERE annotations, 3 levels of biotypes are reported: name, class and family of the ERE overlapping at the genomic location. 

	.. thumbnail:: _images/genomic_ere_annotation.png
      		         
-----------

**full_info_biotypes**
^^^^^^^^^^^^^^^^^^^^^^

.. _Genomic_and_ERE_Annotations_Full:

`1_Genomic_and_ERE_Annotations_Full.csv`: 
For a given peptide carrying an MCS mapped to a given genomic location, BamQuery reports all overlapping biotypes in the Ensembl and Repeat Masker annotations.

Overlap information:

	(a) gene 
	(b) transcript
	(c) genomic location
	(d) ERE name
	(e) ERE class
	(f) ERE family 
	(g) also, the total count of RNA-seq reads bearing the given MCS at the given location in each BAM/CRAM included in **BAM_directories.tsv**.

.. thumbnail:: _images/genomic_and_ERE_Annotations_Full.jpg

-----------

.. _Genomic_and_ERE_Annotations_Summary_Full:

`2_Genomic_and_ERE_Annotations_Summary_Full.csv`: 
For a given peptide at a given genomic location, BamQuery reports all overlapping biotypes in the Ensembl and Repeat Masker annotations.

Overlap information:

	(a) gene 
	(b) transcript
	(c) genomic location
	(d) ERE name
	(e) ERE class
	(f) ERE family 
	(g) also, the total count of RNA-seq reads at the given location in each BAM/CRAM included in **BAM_directories.tsv**.

.. thumbnail:: _images/genomic_and_ERE_Annotations_Summary_Full.jpg

-----------

.. _Genomic_and_ERE_Anno_by_Region_Full:

`3_Genomic_and_ERE_Anno_by_Region_Full.csv`: 

For a given peptide at a given location, BamQuery reports a consensus biotype according to all the overlapping biotypes at the given genomic location. |br|
Each biotype, whether from Ensembl or Repeat Masker annotations, has equal weight in the calculation of the consensus, which is based on the frequencies of the biotype at the location. 

For example: 

One location was collected for a given peptide.
	(a) In the same location the peptide overlaps in frame two transcripts of a canonical protein (``100% in_frame``). 
	(b) In the same location the peptide overlaps the 3'UTR of one transcript of the same canonical protein (``100% 3'UTR``). 
	(c) In the same location the peptide overlaps an ERE region (``100% ERE``).

The final biotype for the peptide at the given location corresponds to: |br|
``In_frame : 50%, 3'UTR : 25%, ERE : 25%``. i.e, Computaiton doesn't take into consideration the transcription expression!.

Best guess : 
			1. 'In_Frame' if it is among the genomic position biotypes.
			2. Otherwise, the biotype with the highest percentage representation in annotation frequencies. If BamQuery detects that all biotypes have equal representation, it will report all of them as the "Best guess".


.. thumbnail:: _images/genomic_and_ERE_Anno_by_Region_Full.jpg

-----------

**summary_info_biotypes**
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _General_Gen_and_ERE_Biotype_Consensus:

`1_General_Gen_and_ERE_Biotype_Consensus.csv`: 
For a given peptide, BamQuery reports assigns a consensus biotype taking into account all possible locations where the peptide appears in the genome. 
To assign the consensus biotype to the peptide, BamQuery aggregates all reads assigned to each biotype that were distributed according to the percentage of the biotype at each genomic location. This distribution is described in more detail in the file 3_Genomic_and_ERE_Anno_by_Region_Full.csv. 
Next, the total reads assigned to each biotype are weighted with the total reads of all samples investigated.

For example: 

Three locations were collected for a given peptide.
	(a) At location 1 the peptide overlaps in-frame with two transcripts of a canonical protein (``100% in_frame``) and at that location there is a count of 10 RNA-seq reads. 
	(b) At location 2 the peptide overlaps in-frame with one transcript and in the 3'UTR of another transcript (``50% in_frame, 50% 3'UTR``) and at that location there are 20 RNA-seq reads. 
	(c) At location 3 the peptide overlaps with the intronic region of a transcript (``100% Intron``) and at that location there are 7 RNA-seq reads.

The final biotype for the peptide corresponds to: |br|
``In_frame : 54%, 3'UTR : 27%, Introns : 19%``. 

Best guess : 
			1. 'In_Frame' if it is among the genomic position biotypes.
			2. Otherwise, the biotype with the highest percentage representation in annotation frequencies. If BamQuery detects that all biotypes have equal representation, it will report all of them as the "Best guess".


.. thumbnail:: _images/general_Gen_and_ERE_Biotype_Consensus.jpg

-----------

.. _Weighted_Gen_and_ERE_Biotype_Consensus:

`2_Weighted_Gen_and_ERE_Biotype_Consensus.csv`: 

For a given peptide, BamQuery reports the consensus of biotypes according to their overlap at all genomic locations of the peptide. 
The biotype representation (percentage) is computed from the count of RNA-seq reads attributed to each biotype according to the coefficients estimated using the ``EM`` algorithm as a function of the total reads for the given peptide in all the samples, (only expressed locations are taken into account to calculate the percentage) .

For example: 

Three locations were collected for a given peptide.
	(a) At location 1 the peptide overlaps in-frame with two transcripts of a canonical protein (``100% in_frame``) and at that location there is a count of 10 RNA-seq reads. 
	(b) At location 2 the peptide overlaps in-frame with one transcript and in the 3'UTR of another transcript (``89% in_frame, 50% 3'UTR``) and at that location there are 20 RNA-seq reads. 
	(c) At location 3 the peptide overlaps with the intronic region of a transcript (``100% Intron``) and at that location there are 7 RNA-seq reads.

The final biotype for the peptide corresponds to: |br|
``In_frame : 75%, 3'UTR : 6%, Introns : 19%``. 

Best guess : 
			1. 'In_Frame' if it is among the genomic position biotypes.
			2. Otherwise, the biotype with the highest percentage representation in the weighted biotype. If BamQuery detects that all biotypes have equal representation, it will report all of them as the "Best guess".

.. thumbnail:: _images/Weighted_Gen_and_ERE_Biotype_Consensus.jpg

-----------

.. _Group_Samples_Gen_and_ERE_Biotype_Consensus:

`3_Group_Samples_Gen_and_ERE_Biotype_Consensus.csv`: 

For a given peptide, BamQuery reports the consensus of biotypes according to their overlap at all genomic locations of the peptide. 
The biotype representation (percentage) is computed from the count of RNA-seq reads attributed to each biotype according to the coefficients estimated using the ``EM`` algorithm as a function of the total reads for the given peptide in every group of samples as well as for all the samples, (only expressed locations are taken into account to calculate the percentage) .

Best guess : 
			1. 'In_Frame' if it is among the genomic position biotypes.
			2. Otherwise, the biotype with the highest percentage representation in the weighted biotype. If BamQuery detects that all biotypes have equal representation, it will report all of them as the "Best guess".

.. thumbnail:: _images/group_Samples_Gen_and_ERE_Biotype_Consensus.jpg

-----------

.. _normal_mode_count_norm_info:

**normal_mode_count_norm_info.xlsx**
------------------------------------

`Sheet : Alignments Read count RNA-seq`

This sheet reports for each peptide queried, all positions in the genome that are perfect alignments for one or more coding sequences of a peptide are reported. For each position, the strand, coding sequence and read count for each BAM/CRAM file are reported.

.. thumbnail:: _images/alignments_Read_count_RNA_seq.jpg

.. _read count RNA seq by peptide:


`Sheet : Read count RNA-seq by peptide`

This sheet reports for each peptide queried, the total reads for each BAM/CRAM file considering all positions. 

.. thumbnail:: _images/read_count_RNA_seq_by_peptide.jpg

.. _log10 RPHM RNA seq by peptide:


`Sheet : log10(RPHM) RNA-seq by peptide`

This sheet reports for each peptide queried, the :math:`rphm` `(read per hundred million)` for each BAM/CRAM file considering all expressed positions. The :math:`rphm = (read\_overlap * 10^8)/total\_primary\_reads` with `total_primary_reads` representing the total number of reads sequenced in a given RNA-Seq experiment. These values are transformed into logarithm :math:`log_{10}(rphm + 1)`.

This information is used to plot the rphm heat map. See `heat maps folder`_

.. thumbnail:: _images/log10RPHM.jpg


======================

.. include:: light_mode.rst

======================

.. _single_cell_example:

.. include:: single_cell.rst

======================

.. include:: translation.rst

.. |br| raw:: html

      <br>