# README for Taff et al.

This repository includes the full set of code and data required to reproduce the analyses, figures, and manuscript of Taff et al. 2023: Joint effects of social interactions and environmental challenges on physiology, internal microbiome, and reproductive performance in tree swallows (*Tachycineta bicolor*). 

Descriptions of files are included below and code scripts also include extensive in line annotation. The repository is permanently archived at Zenodo.



# CODE

There are three scripts that sequentially produce all the outputs.

1. 0_manuscript.Rmd

  - This markdown file produces the formatted manuscript and supplemental material for the paper. Code chunks are included to pull in figures and tables from the other scripts, but only minimal processing is included in the file and it is mostly text. 
  
2. 1_main_r_script.R

 - This script completes most of the analyses presented in the paper. It relies on a number of input data files that are included in the repository and it produces a variety of outputs for tables and figures that are saved to file to be loaded in the manuscript. Note that some tables in the supplementary material are produced by this script but are saved manually as text files to be loaded and formatted in the manuscript.
 
3. 2_microbiome_processing.R 

 - This script must be run first to process the raw sequence files from microbiome samples. In addition to data files included in the repository, running this script will require accessing the raw sequence data from NCBI SRA and downloading the Silva rRNA reference database. The processed output from this script is used as input in the main analysis script.
 
# RAW DATA

The repository includes a variety of raw data inputs that are used for various analyses.

1. 16s_sequences: These sequences are not included in the GitHub repository but are instead archived on NCBI SRA. To fully run the scripts above these read files should be downloaded and added to the repository files.

2. '1_raw_data' folder.

This folder includes all the tabled raw data used in the manuscript in six files.

- *daily_provision_data.txt*: Includes the number of provisioning trips for each day of study for each male and female included in the population. Note that each individual female/nest will have many rows of data per year for different days.

- *data_by_female.txt*: Includes one row for each unique female by year combination with columns for each female measure used in analyses.

- *data_by_micro_sample.txt*: Includes one row for each microbiome swab that was processed with metadata about the bird associated to that swab. Note that each bird will may have multiple rows of microbiome swab data.

- *data_by_nestling.txt*: Includes one row for each nestling included in the study with columns for all nestling measures used.

- *nest_visit_data.txt*: Includes one row for every day of observation with social interaction data from the RFID network. Note that each nest includes multiple days (rows) of data.

- *total_rfid_reads.txt*: Includes the total number of reads (of any bird) at each rfid unit for every day in the study. Note that each nest has many rows of data in each year.

# INTERMEDIATE DATA

In addition to the raw data described above, the repository contains a variety of derived files that are produced by the main analysis script. These include intermediate cleaned data files and saved model or figure objects. For convenience, these files are saved in the repository in the folders '2_modified_data' (data tables) and '3_r_scripts' (figures). They can also be reproduced and modified by running the microbiome processing script and main script again.