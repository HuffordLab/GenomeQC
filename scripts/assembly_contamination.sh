#!/bin/bash

# This script calculates, plots and then email the busco scores for the user uploaded genome assembly and also compares the busco scores with the user selected reference genomes.
# The last part of this script screens the uploaded genome assembly for vector/adaptor sequence contamination.

# Usage: bash assembly_contamination.sh <genome assembly file (fasta format)> <busco output directory name> <busco dataset> <augustus species> <email address> <name for output plots> <reference genome names>


#This command calls the busco script according to the input parameters provided by the user
python /home/nancy/assembly_app/busco-master/scripts/run_BUSCO_fasta.py --in $1 --out $2 -m genome --lineage $3 $4 -c 4 --blast_single_core 
##This command calls a python script which generates and emails a stack plot comparing the scores for the user uploaded genome assembly with the user selected reference genome assembly  
python /home/nancy/assembly_app/busco-master/scripts/fastaplot_busco.py $5 $6 $2 $7
#This command calls a python script that blasts the input genome assembly to the NCBI Univec database and generates a plot showing the number of scaffolds/contigs with blast hits
python /home/nancy/assembly_app/busco-master/scripts/contamination.py $1 $2 $5


