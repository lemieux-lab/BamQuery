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
<h3>3. Install python 3 and create a virtual environment</h3>

<h4>Option 1: Installation with Conda </h4>

For users having no administrator priviledges, we recommend installing BamQuery with conda (link: https://docs.conda.io/en/latest/miniconda.html).

First create a conda environment, activate it, and install all non-python dependencies:

        conda create -n BQ python=3
        conda activate BQ
        conda install -y -c conda-forge r-ggplot2=3.3.6
        conda install -y -c conda-forge r-data.table
        conda install -y -c bioconda bedtools
        conda install -y -c bioconda star=2.7.9a
        
Then create a python virtual environment, activate it, and install all python dependencies:

        cd ${INSTALLDIR}
        python3 -m venv env
        source ${$INSTALLDIR}/env/bin/activate
        pip install --upgrade pip
        pip install pandas
        pip install pysam
        pip install pathos
        pip install xlsxwriter
        pip install seaborn
        pip install billiard
        pip install scipy
        
Launch the analysis:

        conda activate BQ
        source ${$INSTALLDIR}/env/bin/activate
        python3 $INSTALLDIR/BamQuery/BamQuery.py path_to_input_folder name_exp
        


<h4>Option 2: Installation from source</h4>

Python: https://www.python.org/

        python3 -m venv bamquery-venv
        source $INSTALLDIR/bamquery-venv/bin/activate
        
Install python packages in the virtual environment</h4>

        pip install --upgrade pip
        pip install pandas
        pip install pysam
        pip install pathos
        pip install xlsxwriter
        pip install seaborn
        pip install billiard
        pip install numpy
        pip install scipy
        
        
---------------------

**4. Install external dependencies so that their binaries are available in your $PATH:**

STAR 2.7.9a: https://github.com/alexdobin/STAR

bedtools: https://bedtools.readthedocs.io/en/latest/

R: https://www.r-project.org/

Required R packages: ggplot2, data.table

--------------------------

**5. Launch the analysis**

        python3 $INSTALLDIR/BamQuery/BamQuery.py path_to_input_folder name_exp
        




--------------------------


<h2>  Installation using the provided docker container </h2>


A docker container is also available to provide a self contained working environment.

<h3>1. Create an install folder:</h3>

        export INSTALLDIR=/opt/bamquery
        mkdir $INSTALLDIR
        cd $INSTALLDIR

------------------------  

<h3>2. Download the docker image:</h3>

        wget https://bamquery.iric.ca/download/bamquery-2022-12-14.tar.gz
        
------------------------       
<h3>3. Install the docker image (requires sudo access):</h3>
        
        gunzip bamquery-2022-12-14.tar.gz
        sudo docker load --input bamquery-2022-12-14.tar
        
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
        iric/bamquery:0.2 python3 /opt/bamquery/BamQuery/BamQuery.py path_to_input_folder name_exp
        
        
making sure to map any required folder mentionned in the input files (BAM locations, input folder) so that these paths may be available from within the container. 
This is done with multiple arguments `-v $DATAFOLDER:$DATAFOLDER` (where `$DATAFOLDER` is to be replaced by an actual folder name) and `-v $PWD:$PWD` if needed.
Note also that we force the application to run with user permissions instead of root using the `--user $(id -u):$(id -g)` argument.

For more information on configuration, see : https://bamquery.iric.ca/documentation/installation.html
