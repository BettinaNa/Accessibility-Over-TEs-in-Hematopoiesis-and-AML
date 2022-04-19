#!/usr/bin/env bash
#----------------------------------------------------------------
## Run LOLA ReMap 2018 data base and HOMER motif analysis on all relevant subfamilies
## The HOMER run has been removed from this workbook in order to save time. However, using the commands in this workbook the run can be re-generated.
#----------------------------------------------------------------
PATH=$PATH/:/homer/bin/
loladb="/data/lola_data/"
dbname="remap_2018"

lola_runs="LT.HSC
Primitive.only
LSCpos.only
Primitive.LSCpos"

for run in $lola_runs; 
do
    queryfiles=$(less /data/lola_data/$run'_TEs.txt' )

    # create combined_lola file
    mkdir -p /results/LOLA/$run
    echo $'antibody\tqValue\toddsRatio\tsupport\tTE_subfamily' > /results/LOLA/$run/combined_lola.txt

    for query in $queryfiles;
    do
        
        queryname=$(echo $query | awk -F"/" '{print $NF}' | sed -e 's/.bed.sorted//g')
        opname=$queryname"_"$dbname

        mkdir -p /results/LOLA/$run/$opname
        
        # intersect TE file with consensus file
        intersectBed -a $query -b /data/lola_data/Consensus.catalogue.$run.bed -u > /results/LOLA/$run/$opname/querytmp.bed

        # run LOLA
        Rscript /code/helper_scripts/lola.R \
        /results/LOLA/$run/$opname/querytmp.bed \
        /data/lola_data/ALL.TEs.hg38.Overlapping.HEMAT.LSCs.SM.Consensus.bed \
        $loladb \
        /results/LOLA/$run/$opname

        # merge LOLA output from all querynames, qValue<0.05
        awk_script='NR>1 && $22<0.05 {print $18"\t"$22"\t"$5"\t"$6"\t""'$queryname'"}'
        awk "$awk_script" /results/LOLA/$run/$opname/allEnrichments.tsv > /results/LOLA/$run/$opname/$queryname'_allEnrichments_sig.txt'
        cat /results/LOLA/$run/combined_lola.txt /results/LOLA/$run/$opname/$queryname'_allEnrichments_sig.txt' > /results/LOLA/$run/$opname/$queryname'_allEnrichments_sig_comb.txt'
        mv /results/LOLA/$run/$opname/$queryname'_allEnrichments_sig_comb.txt' /results/LOLA/$run/combined_lola.txt

        # # run HOMER
        # the actual homer run has been removed from this workbook, as this is very time consuming. The results can be produced using the folling commands:
        # opname1=$queryname"_MOTIF"
        # mkdir -p /results/HOMER/$run/$opname1
        # /homer/bin/findMotifsGenome.pl /results/LOLA/$run/$opname/querytmp.bed hg38 /results/HOMER/$run/$opname1/ -bg /data/lola_data/ALL.TEs.hg38.Overlapping.HEMAT.LSCs.SM.Consensus.bed -size 200

        rm /results/LOLA/$run/$opname/querytmp.bed
    done

done
