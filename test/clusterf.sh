#!/bin/bash

P=$1

#INPUT=../../vsearch-data/Rfam_9_1.fasta
#INPUT=../../vsearch-data/AF091148.fsa
INPUT=../../vsearch-data/BioMarKs50k.fsa

THREADS=0

USEARCH=$(which usearch)
VSEARCH=../bin/vsearch

if [ "$P" == "u" ]; then
    PROG=$USEARCH
else
    if [ "$P" == "v" ]; then
        PROG=$VSEARCH
    else
        echo You must specify u or v as first argument
        exit
    fi
fi

CMD="/usr/bin/time $PROG \
    --cluster_fast $INPUT \
    --threads $THREADS \
    --id 0.9 \
    --maxaccepts 1 \
    --maxrejects 8 \
    --sizeout \
    --centroids f.$P.centroids \
    --uc f.$P.uc \
    --alnout f.$P.alnout \
    --blast6out f.$P.bl6 \
    --matched f.$P.matched \
    --notmatched f.$P.notmatched \
    --fastapairs f.$P.fastapairs"

#    --sizein
#    --output_no_hits
#    --uc_allhits
#    --query_cov 0.5
#    --clusters files/s.$P.clusters
#    --consout s.$P.consout
#    --msaout s.$P.msaout

echo Cluster test
echo
echo Running command: $CMD
echo

$CMD
