
.. _light_mode_example:

*******************
light mode example
*******************

----------

The light mode of BamQuery was designed to perform a quick search for peptide expression in the specified BAM/CRAM files. In this mode, BamQuery only reports peptide counts and normalization, so no biotyping analysis is performed for peptides (heat map, pie chart are not produced).

**Input**
#########


**Command line:**

.. code::

	BamQuery.py ./normal_mode_example/Input normal_mode_example v38_104 --light


As for normal mode, the input folder `path_to_input_folder` must containt the files : **BAM_directories.tsv** and **peptides.tsv**.

**Output**
##########

BamQuery creates an **output** directory in the same path as the input folder.

This directory contains 3 folders and the main results are organized as follows:

.. code::

	├── alignments
	│   ├── light_mode_info_alignments.xlsx
	│   └── missed_peptides.info
	├── logs
	│  ├── BamQuery_Res_light_mode.log
	│  └── Information_BAM_directories.log
	└── res_light
	    └── light_mode_count_norm_info.xlsx

The output files in `light mode` are similar to those in `normal mode`, see `output_normal_mode_example`_ for detailed information about the output files.

.. note::
		After running BamQuery in `--light` mode, it is possible to run BamQuery in normal mode (to obtain biotype classification and other output files) for a subset of peptides (peptides of interest `PoIs`). To do this, first run BamQuery in light mode by adding the `--dev` parameter. Once BamQuery light has finished, modify the **peptides.tsv** file, to remove the peptides you are no longer interested in. Finally, run the BamQuery search in `normal mode` by removing the `--light` option from the command line. By doing this, BamQuery takes the information already obtained for the expression in light mode and produces heatmap plots and does the biotype analysis only for the `PoIs`.

   .. warning::
   		WARNING: do not modify the **BAM_directories.tsv**, otherwise you will not have consistent information.

