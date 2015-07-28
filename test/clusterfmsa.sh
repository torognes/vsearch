#!/bin/bash

P=$1

#INPUT=../data/Rfam_9_1.fasta
#INPUT=../data/AF091148.fsa
INPUT=../data/BioMarKs50k.fsa

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
    --cluster_fast $INPUT \
    --threads $THREADS \
    --id 0.97 \
    --sizeout \
    --consout s.$P.consout \
    --profile s.$P.profile \
    --msaout s.$P.msaout"

echo Cluster test with msa
echo
echo Running command: $CMD
echo

$CMD
