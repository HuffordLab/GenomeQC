#!/bin/bash

awk '/^>/{print ">" ++i; next}{print}' $1  > $1_renamed
python3 /opt/NG.py $1 $2 $3
python3 /opt/assembly_stats.py $1 $4 $3
python3 /opt/contamination_check/contamination.py $1 $5 $6
/opt/busco/scripts/run_BUSCO.py -i $1 -l $7 -o $8 -m genome -sp $9 -c ${10} 
/opt/gt-1.5.10-Linux_x86_64-64bit-complete/bin/gt suffixerator -db $1_renamed -indexname $1_renamed -tis -suf -lcp -des -ssp -sds -dna
/opt/gt-1.5.10-Linux_x86_64-64bit-complete/bin/gt ltrharvest -index $1_renamed -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 -motif TGCA -motifmis 1 -similar 85 -vic 10 -seed 20 -seqids yes > $1_renamed.harvest.scn
/opt/LTR_FINDER.x86_64-1.0.5/run_ltrfinder_in_pieces.pl -seq $1_renamed -threads 40 -harvest_out -size 1000000 -time 300
cat $1_renamed.harvest.scn $1_renamed.finder.combine.scn > $1_renamed.rawLTR.scn
/LTR_retriever/LTR_retriever -inharvest $1_renamed.rawLTR.scn -threads ${10} -genome $1_renamed 
