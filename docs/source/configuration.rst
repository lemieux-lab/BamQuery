.. _configuration:

##############
Configuration
##############

Bam_files_info.dic
*******************

BamQuery stores in a Python dictionary information about each BAM/CRAM file queried. 
The key for each bam/cram file is obtained from the path where the file is located. 
Therefore, take precautions to store bam/cram files under informative folder names that serve to differentiate them. 
For example, for all GTEx healthy tissue cram files, the file organization should be as follows: 

.. code::

        GTEx
        ├── adipose_subcutaneous
        │   ├── SRR1333352
        │   ├── SRR1338301
        │   ├── SRR1338627
        │   ├── SRR1339740
        │   ├


With this file structure, a cram path file for GTEx should look like this:  ``/home/GTEX/brain_amygdala/SRR1333352/SRR1333352.cram``.

In this example, BamQuery creates the key ``brain_amygdala_SRR1333352`` to save information related to this sample.

This information is organized in a list as follows: |br| 
0: ``/home/GTEX/brain_amygdala/SRR1333352/SRR1333352.cram`` --> Whole path to bam/cram file |br| 
1: ``80302110`` --> Total Primary Read count in the bam/cram file |br| 
2: ``brain_amygdala`` --> Tissue |br| 
3: ``Brain`` --> Tissue type  |br| 
4: ``no`` --> Shortlist |br| 
5: ``NA`` --> Sequencing  |br| 
6: ``NA`` --> Library |br| 
7: ``User_1`` --> The user that includes the bam/cram file information (first user quering a given bam file) |br| 


The Tissue, Tissue type and Shortlist fields must be provided by the first user who queries the given bam/cram file. This is done only once (see instructions below). |br| 

The sequencing and library fields are guessed directly by BamQuery from the bam/cram file. This is also done once when a user configures BamQuery to query the file taking into account its stradedness. 



Provide details to each Bam file
********************************

Every time a BAM file is queried for the first time, you need to provided some information about the origin of the file. |br| 
This is why the following exception will appear when running BamQuery:

.. py:exception:: fill in the `bam_files_tissues.csv` file with the requested information:

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
Yes or No. This sets the BAM file as part of a selected group of samples within a tissue type to calculate the average level of transcript expression.


Once the file :code:`bam_files_tissues.csv` has been filled, you can relaunch BamQuery.

.. |br| raw:: html

      <br>
