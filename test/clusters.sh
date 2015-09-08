#!/bin/bash

P=$1

INPUT=../data/Rfam_9_1.fasta
#INPUT=../data/AF091148.fsa
#INPUT=../data/BioMarKs50k.fsa

THREADS=0

USEARCH=$(which usearch)
VSEARCH=../src/vsearch

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
    --cluster_smallmem $INPUT \
    --usersort \
    --strand both \
    --id 0.9 \
    --maxaccepts 1 \
    --maxrejects 8 \
    --sizeout \
    --centroids s.$P.centroids \
    --uc s.$P.uc \
    --alnout s.$P.alnout \
    --blast6out s.$P.bl6 \
    --matched s.$P.matched \
    --notmatched s.$P.notmatched \
    --consout s.$P.consout \
    --fastapairs s.$P.fastapairs"

#    --sizein
#    --threads $THREADS
#    --minseqlength 1
#    --output_no_hits
#    --uc_allhits
#    --query_cov 0.5
#    --clusters files/s.$P.clusters
#    --msaout s.$P.msaout

echo Cluster test
echo
echo Running command: $CMD
echo

$CMD
