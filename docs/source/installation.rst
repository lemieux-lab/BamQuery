.. _installation:

#############
Installation
#############

                

A. Installation From source
****************************

1. Clone repository from github
###############################

.. code::

        export INSTALLDIR=/opt/bamquery
        mkdir $INSTALLDIR
        cd $INSTALLDIR
        git clone https://github.com/lemieux-lab/BamQuery.git


.. _Install_required_library:

2. Install required library files within $INSTALLDIR:
#####################################################

.. code::

        wget https://bamquery.iric.ca/download/lib_essentials.tar.gz
        tar vxzf lib_essentials.tar.gz

2.a Installation of genomes
^^^^^^^^^^^^^^^^^^^^^^^^^^^
BamQuery supports three different versions of the human genome (v26_88 / v33_99 / v38_104) and two versions of the mouse genome (GRCm38 and GRCm39, respectively: M24 / M30).

You need to download the human or mouse genome version you wish to use to:

.. code::

        cd lib/genome_versions

And use the command below to download any human genome version : v26_88 or v33_99 or v38_104.

.. code::

        wget https://bamquery.iric.ca/download/genome_SET_VERSION.tar.gz
                
or to download any mouse genome version : m24, m30.

.. code::

        wget https://bamquery.iric.ca/download/genome_mouse_SET_VERSION.tar.gz


Finally, you need to:

.. code::

        tar vxzf GENOME_VERSION.tar.gz

2.b Installation of SNPs
^^^^^^^^^^^^^^^^^^^^^^^^
BamQuery supports three different versions of dbSNPs of the human genome (149/151/155) and two versions of dbSNPs of the mouse genome (snps_GRCm38 and snps_GRCm39, respectively: M24 / M30).

You can download the annotated snps you need to (by default BamQuery does not use snps):

.. code:: 

        cd lib/snps

And use the command below to download any dbSNP corresponding to human genome releases : 149 or 151 or 155.

.. code::

        wget https://bamquery.iric.ca/download/dbsnps_SET_RELEASE.tar.gz
                
or to download any dbSNP corresponding to mouse genome releases : GRCm38 or GRCm39.

.. code::

        wget https://bamquery.iric.ca/download/snps_mouse_SET_RELEASE.tar.gz
        
Finally, you need to:

.. code::

        tar vxzf SNPS_RELEASE.tar.gz


3. Create a virtual environment and install dependencies
#########################################################

Option 1: Installation with Conda
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For users having no administrator priviledges, we recommend installing BamQuery with conda: https://docs.conda.io/en/latest/miniconda.html |br|

1. First create a conda environment and activate it: |br|

.. code::

        conda create -n BQ
        conda activate BQ

2. Then install all dependencies:

2. Then install all dependencies: |br|

.. code::
        
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

        
3. Launch the analysis: |br|

.. code::

        conda activate BQ
        source ${INSTALLDIR}/env/bin/activate
        python3 ${INSTALLDIR}/BamQuery/BamQuery.py path_to_input_folder name_exp genome_version
        

Option 2: Installation from source
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. Download Python 3 and creare a virtual environment. Python: https://www.python.org/ |br|

.. code::

        python3 -m venv bamquery-venv
        source ${INSTALLDIR}/bamquery-venv/bin/activate

2. Install python packages in the virtual environment |br|

.. code::

        pip install --upgrade pip
        pip install pandas
        pip install pysam
        pip install pathos
        pip install xlsxwriter
        pip install seaborn
        pip install billiard
        pip install numpy
        pip install scipy
        

3. Install external dependencies so that their binaries are available in your $PATH:

STAR 2.7.9a: https://github.com/alexdobin/STAR |br|
bedtools: https://bedtools.readthedocs.io/en/latest/ |br|
R: https://www.r-project.org/, required R packages: ggplot2, data.table |br|


4. Launch the analysis

.. code::

        python3 ${INSTALLDIR}/BamQuery/BamQuery.py path_to_input_folder name_exp genome_version



=======================


B. Installation using the provided docker container
***************************************************

A docker container is also available to provide a self contained working environment. |br|

1. Create an install folder:
############################

.. code::

        export INSTALLDIR=/opt/bamquery
        mkdir $INSTALLDIR
        cd $INSTALLDIR

2. Download the docker image:
#############################

.. code::

        wget https://bamquery.iric.ca/download/bamquery-2023-02-17.tar.gz

3. Install the docker image (requires sudo access):
###################################################

.. code::

        gunzip bamquery-2023-02-17.tar.gz
        sudo docker load --input bamquery-2023-02-17

4. Install required library files within $INSTALLDIR:
#####################################################

Please, follow the instructions in step 2 enumerated above. See :ref:`Install_required_library`

5. Launch the analysis from the docker container:
#################################################

.. code::

        sudo docker run -i -t  \
        --user $(id -u):$(id -g) \
        -v $INSTALLDIR/lib:/opt/bamquery/lib \
        -v $DATAFOLDER:$DATAFOLDER  \
        -v $PWD:$PWD \
        iric/bamquery:0.2 python3 /opt/bamquery/BamQuery/BamQuery.py path_to_input_folder name_exp

making sure to map any required folder mentionned in the input files (BAM locations, input folder) so that these paths may be available from within the container.  This is done with multiple arguments :code:`-v $DATAFOLDER:$DATAFOLDER` (where :code:`$DATAFOLDER` is to be replaced by an actual folder name) and :code:`-v $PWD:$PWD` if needed. |br|
Note also that we force the application to run with user permissions instead of root using the :code:`--user $(id -u):$(id -g)` argument.

For more information on configuration, see :ref:`configuration`.

.. note::

        BamQuery requires a specific folder structure to work. |br|
        Once BamQuery is installed, check that the structure looks as follows:

        .. code-block::

                .
                ├── BamQuery
                │   ├── BamQuery.py
                │   ├── genomics
                │   ├── plotting
                │   ├── readers
                │   ├── README.md
                │   └── utils
                └── lib
                    ├── coefficients.dic
                    ├── Cosmic_info.dic
                    ├── ERE_info.dic
                    ├── ERE_info_mouse.dic
                    ├── EREs_souris.bed
                    ├── genome_versions
                    │   ├── genome_mouse_m24
                    │   ├── genome_mouse_m30
                    │   ├── genome_v26_88
                    │   ├── genome_v33_99
                    │   └── genome_v38_104
                    ├── hg38_ucsc_repeatmasker.gtf
                    ├── README.txt
                    └── snps
                        ├── snps_dics_149
                        ├── snps_dics_149_common
                        ├── snps_dics_151
                        ├── snps_dics_151_common
                        ├── snps_dics_155
                        └── snps_dics_155_common

.. |br| raw:: html

      <br>
