#!/bin/bash

INPUT=../data/BioMarKs.fsa
#INPUT=../data/PR2-18S-rRNA-V4.fsa
#INPUT=../data/Rfam_9_1.fasta

THREADS=0

VSEARCH=../src/vsearch
USEARCH=$(which usearch)

OUTDIR=.

MINSIZE=1
MAXSIZE=10000000
TOPN=1000000

P=$1

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


echo Creating test database

rm -rf temp.fsa
for n in 1 2 3 4 5 6 7 8 9 10; do
#for n in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do
    cat $INPUT >> $OUTDIR/temp.fsa
done

$VSEARCH --shuffle $OUTDIR/temp.fsa --output $OUTDIR/temp2.fsa

INPUT=$OUTDIR/temp2.fsa

echo 

CMD="/usr/bin/time $PROG \
        --derep_prefix $INPUT \
        --threads $THREADS \
        --minuniquesize $MINSIZE \
        --topn $TOPN \
        --sizeout \
        --output $OUTDIR/derep.$1.fsa \
        --uc $OUTDIR/derep.$1.uc"

echo
echo Prefix dereplication test
echo
echo Running: $CMD
echo

$CMD
