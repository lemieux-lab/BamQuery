# BamQuery: a proteogenomic tool to explore the immunopeptidome and prioritize actionable tumor antigens

BamQuery is a computational pipeline that counts all RNA-seq reads that can code for a given peptide (8 to 11 residues) in chosen RNA-seq samples (single or bulk). The strength of BamQuery lies in its ability to quantify RNA expression from any genomic region: protein-coding exons or non-coding genomic regions, annotated or not. Briefly, it reverse-translates peptides into all possible coding sequences, aligns these sequences on the genome (human or mouse) with STAR, retains only regions (spliced or not) that have perfect alignments with the coding sequence, and queries (grep) the coding sequences at their respective coding regions in the RNA-seq sample (bam file) of interest. All primary reads (samtools view -F 256) that can code for the peptide are counted, normalized on the total primary read number of the sample, and listed in a detailed report.

BamQuery also supports single nucleotide mutated peptides by optionally including the full dbSNP database in the genomic alignment step. Mutated peptides that derive from indels or from regions not listed in dbSNP can be analyzed with the manual mode of BamQuery by providing the genomic region of origin of the peptide. Finally, BamQuery uses an expectation-maximization algorithm to annotate the biotype of each peptide analyzed based on the annotations of their genomic regions of origin and the number of reads present at each of these regions.

BamQuery can be installed as described here or tested with our [online portal](https://bamquery.iric.ca/search). The portal will quantify any peptide in medullary thymic epithelial cells and dendritic cells to evaluate their probability of being expressed by normal tissues. A prediction of their immunogenicity will also be provided based on these expression measures. Therefore, BamQuery can be used to evaluate the specificity and immunogenicity of tumor antigens.

The preprint of our article can be found on BioRXiv: https://doi.org/10.1101/2022.10.07.510944

For detailed usage instructions, see our website: https://bamquery.iric.ca/

BamQuery was designed and developed by Gregory Ehx (GIGA Institute, University of Liege) and Maria Virginia Ruiz Cuevas (Institute for Research in Immunology and Cancer (IRIC), University of Montreal).

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
* [matplotlib](https://pypi.org/project/matplotlib/)
* [pathos](https://pypi.org/project/pathos/)
* [seaborn](https://seaborn.pydata.org/installing.html)
* [scipy](https://scipy.org/install/)

#### Required R package:
* [ggplot2](https://github.com/tidyverse/ggplot2)
* [data.table](https://github.com/Rdatatable/data.table)


__________________________________________________________________________________________

<h1>  Installation </h1>



<h2>  Installation From source </h2>


See the [user manual](https://bamquery.iric.ca/documentation/) for a detailed description usage.

<h3>1. Clone repository from github</h3>

        export INSTALLDIR=/opt/bamquery
        mkdir $INSTALLDIR
        cd $INSTALLDIR
        git clone https://github.com/lemieux-lab/BamQuery.git
        
------------------

<h3>2. Install required library files within $INSTALLDIR:</h3>

        wget https://bamquery.iric.ca/download/lib_essentials.tar.gz
        tar vxzf lib_essentials.tar.gz

 <br>
 
<h4>2.a Installation of genomes</h4>

BamQuery supports three different versions of the human genome (v26_88 / v33_99 / v38_104) and two versions of the mouse genome (GRCm38 and GRCm39, respectively: M24 / M30).

You need to download the human or mouse genome version you wish to use to:
        
        cd lib/genome_versions

And use the command below to download any human genome version : v26_88 or v33_99 or v38_104.

        wget https://bamquery.iric.ca/download/genome_SET_VERSION.tar.gz
                
or to download any mouse genome version : m24, m30.

        wget https://bamquery.iric.ca/download/genome_mouse_SET_VERSION.tar.gz


Finally, you need to:
        
        tar vxzf GENOME_VERSION.tar.gz
        
 <br>
 
<h4>2.b Installation of SNPs</h4>

BamQuery supports three different versions of dbSNPs of the human genome (149/151/155) and two versions of dbSNPs of the mouse genome (snps_GRCm38 and snps_GRCm39, respectively: M24 / M30).

You can download the annotated snps you need to (by default BamQuery does not use snps):
        
        cd lib/snps

And use the command below to download any dbSNP corresponding to human genome releases : 149 or 151 or 155.

        wget https://bamquery.iric.ca/download/dbsnps_SET_RELEASE.tar.gz
                
or to download any dbSNP corresponding to mouse genome releases : GRCm38 or GRCm39.

        wget https://bamquery.iric.ca/download/snps_mouse_SET_RELEASE.tar.gz
        
Finally, you need to:
        
        tar vxzf SNPS_RELEASE.tar.gz

-----------------------
<h3>3. Create a virtual environment and install dependencies</h3>

<h4>Option 1: Installation with Conda </h4>

For users having no administrator priviledges, we recommend installing BamQuery with [conda](https://docs.conda.io/en/latest/miniconda.html).

1. First create a conda environment and activate it:

        conda create -n BQ
        conda activate BQ
        
2. Then install all dependencies:

        conda install -y -c bioconda pysam
        conda install -y -c anaconda pandas
        conda install -y -c conda-forge pathos
        conda install -y -c conda-forge xlsxwriter
        conda install -y -c anaconda seaborn
        conda install -y -c conda-forge billiard
        conda install -y -c anaconda scipy
        conda install -y -c bioconda bedtools
        conda install -y -c bioconda star=2.7.9a
        conda install -y -c conda-forge mamba
        mamba install -y -c conda-forge r-ggplot2
        mamba install -y -c conda-forge r-data.table
        
3. Launch the analysis:

        conda activate BQ
        python3 ${INSTALLDIR}/BamQuery/BamQuery.py path_to_input_folder name_exp genome_version
        


<h4>Option 2: Installation from source</h4>

1. Download Python 3 and creare a virtual environment. Python: https://www.python.org/

        python3 -m venv bamquery-venv
        source ${INSTALLDIR}/env/bin/activate
        
2. Install python packages in the virtual environment</h4>

        pip install --upgrade pip
        pip install pandas
        pip install pysam
        pip install pathos
        pip install xlsxwriter
        pip install seaborn
        pip install billiard
        pip install numpy
        pip install scipy
        

3. Install external dependencies so that their binaries are available in your $PATH:</h3>

* STAR 2.7.9a: https://github.com/alexdobin/STAR

* bedtools: https://bedtools.readthedocs.io/en/latest/

* R: https://www.r-project.org/, required R packages: ggplot2, data.table


4. Launch the analysis</h3>

        python3 ${INSTALLDIR}/BamQuery/BamQuery.py path_to_input_folder name_exp genome_version
        


--------------------------


<h2>  Installation using the provided docker container </h2>


A docker container is also available to provide a self contained working environment.

<h3>1. Create an install folder:</h3>

        export INSTALLDIR=/opt/bamquery
        mkdir $INSTALLDIR
        cd $INSTALLDIR

------------------------  

<h3>2. Download the docker image:</h3>

        wget https://bamquery.iric.ca/download/bamquery-2022-12-22.tar.gz
        
------------------------       
<h3>3. Install the docker image (requires sudo access):</h3>
        
        gunzip bamquery-2022-12-22.tar.gz
        sudo docker load --input bamquery-2022-12-22.tar
        
------------------------        
<h3>4. Install required library files within $INSTALLDIR:</h3>

Please, follow the instructions in step 2 enumerated above. 


------------------------        
<h3>5. Launch the analysis from the docker container:</h3>

        sudo docker run -i -t  \
        --user $(id -u):$(id -g) \
        -v $INSTALLDIR/lib:/opt/bamquery/lib \
        -v $DATAFOLDER:$DATAFOLDER  \
        -v $PWD:$PWD \
        iric/bamquery:0.2 python3 /opt/bamquery/BamQuery/BamQuery.py path_to_input_folder name_exp genome_version
        
        
making sure to map any required folder mentionned in the input files (BAM locations, input folder) so that these paths may be available from within the container. 
This is done with multiple arguments `-v $DATAFOLDER:$DATAFOLDER` (where `$DATAFOLDER` is to be replaced by an actual folder name) and `-v $PWD:$PWD` if needed.
Note also that we force the application to run with user permissions instead of root using the `--user $(id -u):$(id -g)` argument.

For more information on configuration, see : https://bamquery.iric.ca/documentation/installation.html
