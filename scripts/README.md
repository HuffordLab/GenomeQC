This directory contains python and bash scripts which are called by the main script “server.R” to calculate various assembly and annotation metrics. 

**Brief description of all the scripts in this directory:**

- **NG.py:** 

  This script calculates the NG values at various thresholds (1-100%) for the given input genome assembly.

  Usage: python NG.py \<genome assembly file\> <output file name\> <estimated genome size\>


- **annotation.sh:** 

  This script calculates, plots and then email the busco scores for the user uploaded genome annotations and compare those   scores with the user selected reference genomes.

  Usage: bash annotation.sh \<transcripts file\> <busco output directory name\> <busco dataset\> <email address\> <name for output plot\> <reference genome names\>

- **assembly_contamination.sh:**

  This script calculates, plots and then email the busco scores for the user uploaded genome assembly and also compares the busco scores with the user selected reference genomes.

  The last part of this script screens the uploaded genome assembly for vector/adaptor sequence contamination.

  Usage: bash assembly_contamination.sh \<genome assembly file\> <busco output directory name\> <busco dataset\> <augustus species\> <email address\> <name for output plots\> <reference genome names\>

- **assembly_stats.py:**

  This script calculates all the basic length metrics for the given input genome assembly.

  Usage: python assembly_stats.py \<genome assembly file\> <output file name\> <estimated genome size\> 

- **contamination.py:** 

  This script blasts the input genome assembly to the NCBI Univec database and generates a plot showing the number of scaffolds/contigs with blast hits

  Usage: python contamination.py \<genome assembly file\> <output file name\> <email address\> 

- **fastaplot_busco.py:**

  This script plots and email the scores for the user uploaded genome assembly and compare the scores to the user selected reference genomes

  Usage: python fastaplot_busco.py \<email address\> <name for output plot\> <busco output directory name\> <reference genome names\>

- **gff3_stats.py:**

  This script calculates all the basic metrics for the given input genome annotation set.

  Usage: python gff3_stats.py \<genome annotation file\> <output file name\>  

- **gffplot_busco.py:**

  This script plots and email the scores for the user uploaded genome annotations and compare the scores to the user selected reference genomes

  Usage: python gffplot_busco.py \<email address\> <name for output plot\> <busco output directory name\> <reference genome names\>

- **no_fasta_gff.py:**
 
  This script plots and email the scores for the user selected reference genome assemblies and annotations

  Usage: python no_fasta_gff.py \<email address\> <busco assembly output directory name\> <busco annotation output directory name\> <reference genome names\>

- **no_ref_fastaplot.py:**

  This script plots and email the scores for the user uploaded genome annotation set 

  Usage: python no_ref_gffplot.py \<email address\> <name for output plot\> <busco output directory name\> 

- **no_ref_gffplot.py:**
 
  This script plots and email the scores for the user uploaded genome annotation set 

  Usage: python no_ref_gffplot.py \<email address\> <name for output plot\> <busco output directory name\> 

- **noref_annotation.sh:**

  This script calculates, plots and then email the busco scores for the user uploaded genome annotations.

  Usage: bash annotation.sh \<transcripts file\> <busco output directory name\> <busco dataset\> <email address\> <name for output plots\> 

- **noref_assembly_contamination.sh:**

  This script calculates, plots and then email the busco scores for the user uploaded genome assembly.
  The last part of this script screens the uploaded genome assembly for vector/adaptor sequence contamination.

  Usage: bash assembly_contamination.sh \<genome assembly file\> <busco output directory name\> <busco dataset\> <augustus species\> <email address\> <name for output plots\> 

- **run_BUSCO_fasta.py:**

  This script calculates the BUSCO scores for the uploaded genome assembly. 

  Usage: command for running this script is provided in the assembly_contamination.sh and noref_assembly_contamination.sh scripts.

- **run_BUSCO_gff.py:**
  
  This script calculates the BUSCO scores for the uploaded genome annotations. 

  Usage: command for running this script is provided in the annotation.sh and noref_annotation.sh scripts.

