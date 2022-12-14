.. _instalation:

##############################
Installation and configuration
##############################



A. Installation From source
****************************

1. Clone repository from github
###############################

.. code::

        export INSTALLDIR=/opt/bamquery
        mkdir $INSTALLDIR
        cd $INSTALLDIR
        git clone https://github.com/lemieux-lab/BamQuery.git

2. Install required library files within $INSTALLDIR:
#####################################################

.. code::

        wget https://bamquery.iric.ca/download/bamquery-lib.tar.gz
        tar vxzf bamquery-lib.tar.gz

3. Install python 3 and create a virtual environment
####################################################
Python: https://www.python.org/ |br|

.. code::

        python3 -m venv bamquery-venv
        source $INSTALLDIR/bamquery-venv/bin/activate

3.a. Install python packages in the virtual environment |br|

.. code::

        pip install --upgrade pip
        pip install pandas
        pip install pysam
        pip install pathos
        pip install xlsxwriter
        pip install seaborn
        pip install billiard
        pip install plotnine
        pip install sklearn

4. Install external dependencies so that their binaries are available in your $PATH:
####################################################################################

STAR 2.7.9a: https://github.com/alexdobin/STAR |br|
bedtools: https://bedtools.readthedocs.io/en/latest/ |br|
R: https://www.r-project.org/ |br|

Required R packages: ggplot2, data.tables, ggpubr |br|

5. Launch the analysis
######################

.. code::

        python3 $INSTALLDIR/BamQuery/BamQuery.py path_to_input_folder name_exp


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

        wget https://bamquery.iric.ca/download/bamquery-2022-12-06.tar.gz

3. Install the docker image (requires sudo access):
###################################################

.. code::

        gunzip bamquery-2022-12-06.tar.gz
        sudo docker load --input bamquery-2022-12-06.tar

4. Install required library files within $INSTALLDIR:
#####################################################

.. code::

        wget https://bamquery.iric.ca/download/bamquery-lib.tar.gz
        tar vxzf bamquery-lib.tar.gz

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


=======================


Configuration
*************

Every time a BAM file is going to be queried for the first time, BamQuery is going to need information about its origin. |br| 
This is why the following exception will appear when running BamQuery:

.. py:exception:: fill in the :code:`bam_files_tissues.csv` file with the requested information:

    Before to continue you must provide the tissue type for the bam files annotated in the file : .../output/res/AUX_files/bam_files_tissues.csv. Please enter for each sample : tissue, tissue_type, shortlist.

To resolve this, you must fill in the :code:`bam_files_tissues.csv` file with the requested information. |br| 
BamQuery stores the information, so this is a one-time operation for each BAM file. |br| 

Columns in :code:`bam_files_tissues.csv` : |br| 

For each BAM file, you must provide tissue, tissue_type, shortlist. |br| 
This classification is used by BamQuery for the elaboration of the heatmaps. See :ref:`heat maps folder`

**tissue:**
Refers to the tissue of the sample. For example: prostate

**tissue_type:**
It refers to a specific feauture of the tissue. For example: prostate tissue, can be classified as a type of SexSpecific tissue

**shortlist:**
Yes or No. This sets the BAM file as part of the selected samples to calculate the average level of transcript expression associated with the tissue type.


Once the file :code:`bam_files_tissues.csv` has been filled, you can relaunch BamQuery.

.. |br| raw:: html

      <br>