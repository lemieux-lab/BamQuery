.. BamQuery documentation master file, created by
   sphinx-quickstart on Fri Jun  4 20:31:41 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. image:: _images/favicon.png
   :alt: BamQuery logo


BamQuery: a proteogenomic tool to explore the immunopeptidome and prioritize actionable tumor antigens 
======================================================================================================

.. image:: https://img.shields.io/badge/python-3.6-blue.svg 

BamQuery is on `Github`_ and can be used from its `web interface`_.

.. _Github: https://github.com/lemieux-lab/BamQuery

.. _web interface: https://bamquery.iric.ca/

.. _online portal: https://bamquery.iric.ca/search


A Quick Intro:
-----------------

BamQuery is a computational pipeline that counts all RNA-seq reads that can code for a given peptide (8 to 11 residues) in chosen RNA-seq samples (single or bulk). The strength of BamQuery lies in its ability to quantify RNA expression from any genomic region: protein-coding exons or non-coding genomic regions, annotated or not. Briefly, it reverse-translates peptides into all possible coding sequences, aligns these sequences on the genome (human or mouse) with STAR, retains only regions (spliced or not) that have perfect alignments with the coding sequence, and queries (grep) the coding sequences at their respective coding regions in the RNA-seq sample (bam file) of interest. All primary reads (samtools view -F 256) that can code for the peptide are counted, normalized on the total primary read number of the sample, and listed in a detailed report.

BamQuery supports only genomically templated peptides, not proteasome-recombined peptides. Any non-mutated peptide (located on splicing sites or not) is supported. BamQuery supports single nucleotide mutated peptides by optionally including the full dbSNP database in the genomic alignment step. Mutated peptides that derive from indels or from regions not listed in dbSNP can be analyzed with the manual mode of BamQuery by providing the genomic region of origin of the peptide. Finally, BamQuery uses an expectation-maximization algorithm to annotate the most-likely biotype of each peptide analyzed based on the annotations of their genomic regions of origin and the number of reads present at each of these regions.

BamQuery can be installed as described here or tested with our `online portal`_. The portal will quantify any peptide in medullary thymic epithelial cells and dendritic cells to evaluate their probability of being expressed by normal tissues. A prediction of their immunogenicity will also be provided based on these expression measures. Therefore, BamQuery can be used to evaluate the specificity and immunogenicity of tumor antigens.

The preprint of our article can be found on BioRXiv: https://doi.org/10.1101/2022.10.07.510944

BamQuery was designed and developed by Gregory Ehx (GIGA Institute, University of Liege) and Maria Virginia Ruiz Cuevas (*Institute for Research in Immunology and Cancer*  (IRIC_), University of Montreal).


.. _Maria Virginia Ruiz Cuevas: https://github.com/VirginieR
.. _IRIC: http://www.iric.ca

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   method.rst
   command.rst
   examples.rst
   biotype_classification.rst
   installation.rst


