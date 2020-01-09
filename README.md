# GenomeQC: Genome Assembly and Annotation Metrics

GenomeQC generates descriptive summaries with intuitive graphics for
genome assembly and structural annotations. It also benchmarks user
supplied assemblies and annotations against the publicly available
reference genomes of their choice. It is optimized for small and medium
sized genomes (&lt;2.5 Gb) and has pre-computed results for
several maize genomes.

**There is a Dockerfile available (with the associated scripts) to run the pipeline without installing any dependencies.**


**Installation**

**Bioinformatics software dependencies**

**GenomeQC web application calls upon the following bioinformatics tools and
database to perform computation. These tools needs to be installed and
configured in the path of the working directory.**

**At the time of release, this application was tested with:**

-  BUSCO v3.0.2
([*https://gitlab.com/ezlab/busco*](https://gitlab.com/ezlab/busco))
software and its dependencies

-   NCBI BLAST+ v2.28.0
     ([*https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+*](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+)),

-   HMMER v3.1b2 (*<http://hmmer.org/>)*

-   Augustus v3.2.1
     ([*http://bioinf.uni-greifswald.de/augustus/*](http://bioinf.uni-greifswald.de/augustus/))


-   Gffread 0.9.12
     (http://ccb.jhu.edu/software/stringtie/gff.shtml\#gffread\_dl)

-   NCBI UniVec Database (ftp://ftp.ncbi.nlm.nih.gov/pub/UniVec/)
-   Blobtools v1.1 taxify (https://blobtools.readme.io/docs/taxify)

**GenomeQC components:**

GenomeQC is a collection of R and Python scripts. These R scripts need
to be placed in the directory of R Shiny package.

The two main scripts necessary to run the application are ui.R and
server.R.

**ui.R :** This script is the source of user interface definition which
lays out the user interface.

**server.R:** This script, which can be found in the scripts folder of
the GenomeQC Github repository, calls various packages and python and
bash scripts for calculating different metrics.

Running GenomeQC requires a Linux server, R shiny (version 1.5.9) and
Python (version 3.6). Furthermore, it requires the following packages:

  
  
 | R packages     |               |               |                                   
 | -------------  | ------------  |  ------------ |
 | tools          |  R.utils      |  shinyWidgets |   DT
 | seqinr         |  tidyverse    |  shinyBS      |  promises
 | Biostrings     |  gridExtra    |  reshape      |  future
 | stringr        |  grid         |  cowplot      |  

 | Python packages       |             |                          |
 | --------------------- | ------------| ------------------------ | 
 | sys                   |  traceback  |   Bio.Blast.Applications |  email.mime.text
 | os                    | subprocess  |  iglob                   | email.mime.application
 | Bio                   | Statistics  |  pandas                  | email.mime.multipart
 | re                    | Numpy       |  plotly.offline          | smtplib
 | argparse              | collections |  plotly.graph\_objs      | 

**Operating Instructions**

Three modes are available:

**Compare reference genomes:**

 This section outputs various pre-computed assembly and annotation
 metrics from a user-selected list of reference genomes.

**Analyze your genome assembly:**

 This section provides the user the option to perform analysis on their
 genome assembly as well as benchmark their analysis with the
 pre-computed reference genomes.

**Analyze your genome annotations:**

 This section provides the user the option to perform analysis on their
 genome annotations as well as benchmark their analysis with the
 pre-computed reference annotations.
 

**See also an online version of the manual for more details:**
GenomeQC\_userguide.pdf

**Licensing**

GNU GPL V3.

**How to cite GenomeQC**

GenomeQC: A quality assessment tool for genome assemblies and gene structure annotations
Nancy Manchanda, John L. Portwood II, Margaret R. Woodhouse, Arun S. Seetharam, Carolyn J. Lawrence-Dill, Carson M. Andorf, Matthew B. Hufford

bioRxiv 795237; doi: https://doi.org/10.1101/795237

**Acknowledgements**

Funding: This work was supported by the United State Department of
Agriculture (USDA).

**Please send questions to: john.portwood@ars.usda.gov**
