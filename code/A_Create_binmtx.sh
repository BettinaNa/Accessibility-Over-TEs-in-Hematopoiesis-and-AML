#!/usr/bin/env bash
# make binary matrices for all runs
consensusfiles=$( ls /data/consensus_files )

for f in $consensusfiles ;
do
   ##for sample files
   Rscript /code/helper_scripts/createbinarymat.R \
   /data/consensus_files/$f \
   ".narrowPeak" \
   "$f".Binarymat

   ## For TEs 
   Rscript /code/helper_scripts/createbinarymat.R \
   /data/consensus_files/$f \
   ".bed.sorted" \
   "$f".Binarymat.repeats \
   /data/TE_files/ 

done
