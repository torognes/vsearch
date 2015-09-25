#!/bin/bash

P=$1

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
    --cluster_otus $INPUT \
    --strand both \
    --id 0.97 \
    --maxaccepts 1 \
    --maxrejects 8 \
    --sizeout \
    --centroids o.$P.centroids \
    --uc o.$P.uc \
    --alnout o.$P.alnout \
    --blast6out o.$P.bl6 \
    --matched o.$P.matched \
    --notmatched o.$P.notmatched \
    --fastapairs o.$P.fastapairs \
    --otus o.$P.otus"

#    --sizein
#    --threads $THREADS
#    --minseqlength 1
#    --output_no_hits
#    --uc_allhits
#    --query_cov 0.5
#    --clusters files/s.$P.clusters
#    --consout s.$P.consout
#    --msaout s.$P.msaout

echo Testing cluster_otus
echo
echo Running command: $CMD
echo

$CMD
