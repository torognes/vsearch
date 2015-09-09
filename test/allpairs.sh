#!/bin/bash

P=$1

DB=../data/Rfam_11_0.repr.fasta

T=0
ID=0.7

USEARCH=$(which usearch)
VSEARCH=../bin/vsearch

if [ "$P" == "u" ]; then
    PROG=$USEARCH
    THR=""
else
    if [ "$P" == "v" ]; then
        PROG=$VSEARCH
        THR="--threads $T"
    else
        echo You must specify u or v as first argument
        exit
    fi
fi



/usr/bin/time $PROG \
    --allpairs_global $DB \
    $THR \
    --id $ID \
    --query_cov 0.5 \
    --uc a.$P.uc \
    --userfields query+target \
    --userout a.$P.userout.txt \
    --alnout a.$P.alnout.txt \
    --fastapairs a.$P.fastapairs.fsa \
    --blast6out a.$P.blast6.txt \
    --matched a.$P.matched.fsa \
    --notmatched a.$P.notmatched.fsa
