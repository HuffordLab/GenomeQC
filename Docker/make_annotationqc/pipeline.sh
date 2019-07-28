#!/bin/bash

python3 /opt/gff3_stats.py $1 $2 
/opt/busco/scripts/run_BUSCO.py --in $3 --out $4 -m transcriptome --lineage $5 -c $6
 
 
