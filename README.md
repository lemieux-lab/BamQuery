# BamQuery: a proteogenomic tool for the genome-wide exploration of the immunopeptidome


MHC class I–associated peptides (MAPs), collectively referred to as the immunopeptidome, define the immune self for CD8+ T cells and have a pivotal role in cancer immunosurveillance. While MAPs were long thought to be solely generated by the degradation of canonical proteins, recent advances in the field of proteogenomics (genomically-informed proteomics) have evidenced that ∼10% of MAPs originate from allegedly noncoding genomic sequences. Among these sequences, the endogenous retroelements (EREs) are notably under intense scrutiny as possible cancer-specific antigen (TSAs) source.

With the increasing number of cancer-oriented immunopeptidomic and proteogenomic studies comes the need to accurately attribute an RNA expression level to each MAP identified by mass-spectrometry. Here, we introduce BamQuery (BQ), a computational tool to count all reads able to code for any MAP in any RNA-seq data chosen by the user as well as to annotate each MAP with all available biological features. Using BQ, we found that most canonical MAPs can derive from an average of two different genomic regions, whereas most tested ERE-derived MAPs can be generated by numerous (median of 682) different genomic regions and RNA transcripts.

We show that published ERE MAPs considered as TSAs candidates can be coded by numerous other genomic regions than those previously studied, resulting in high undetected expression in normal tissues. Similarly, we also show that some mutated neoantigens previously published as presumably specific anti-cancer targets can in fact be generated by other non-mutated, non-coding, widely expressed RNA-seq reads in normal tissues. In light of these observations, we conclude that BQ could become an essential tool in any TSA-identification/validation pipelines in the near future.

BamQuery is developed by Maria Virginia Ruiz Cuevas at the Institute for Research in Immunology and Cancer (IRIC).

#### Required software:
* [Python 3.6.8](https://www.python.org/downloads/release/python-368/)
* [R 4.2.0](https://cran.r-project.org/src/base/R-4/R-4.2.0.tar.gz)
* [BEDtools](https://bedtools.readthedocs.io/en/latest/)
* [STAR](https://github.com/alexdobin/STAR/archive/2.7.9a.tar.gz)


#### Required Python package:
* [pandas](https://pypi.org/project/pandas/)
* [numpy](https://pypi.org/project/numpy/)
* [billiard](https://pypi.org/project/billiard/)
* [pysam](https://pypi.org/project/pysam/)
* [multiprocessing](https://pypi.org/project/multiprocessing/)
* [subprocess](https://pypi.org/project/subprocess/)
* [math](https://pypi.org/project/math/)
* [matplotlib](https://pypi.org/project/matplotlib/)
* [collections](https://pypi.org/project/collections/)
* [pathos](https://pypi.org/project/pathos/)
* [seaborn](https://seaborn.pydata.org/installing.html)
* [scipy](https://scipy.org/install/)

#### Required R package:
* [ggplot2](https://github.com/tidyverse/ggplot2)
* [data.table](https://github.com/Rdatatable/data.table)

## Installation

### Install via Docker
Docker image of BamQuery is at .
See the [user manual](https://bamquery.iric.ca/documentation/) for a detailed description usage.

### Install from source
1. Install all software listed above.

2. Download or clone the BamQuery repository to your local system:

        git clone (https://github.com/lemieux-lab/BamQuery)

3. Obtain the reference files.

