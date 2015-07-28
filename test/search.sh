#!/bin/bash

P=$1

Q=../data/Rfam_9_1.fasta
DB=../data/Rfam_9_1.fasta
#Q=../data/BioMarKs50k.fsa
#DB=../data/BioMarKs50k.fsa
T=0
ID=0.5
MA=1
MR=32

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
    --usearch_global $Q \
    --db $DB \
    --threads $T \
    --strand plus \
    --id $ID \
    --self \
    --wordlength 8 \
    --sizeout \
    --output_no_hits \
    --gapopen 20I/2E \
    --gapext 2I/1E \
    --maxaccepts $MA \
    --maxrejects $MR \
    --alnout alnout.$P.txt \
    --fastapairs fastapairs.$P.fsa \
    --dbmatched dbmatched.$P.fsa \
    --dbnotmatched dbnotmatched.$P.fsa \
    --blast6out blast6out.$P.bl6"

echo Search test
echo
echo Running command: $CMD
echo

$CMD
