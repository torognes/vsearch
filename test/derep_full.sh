#!/bin/bash

INPUT=../../vsearch-data/PR2-18S-rRNA-V4.fsa
#INPUT=../../vsearch-data/Rfam_9_1.fasta

THREADS=0

VSEARCH=../bin/vsearch
USEARCH=$(which usearch)

OUTDIR=.

MINSIZE=1
MAXSIZE=100000
TOPN=10000

echo Creating test database

rm -rf temp.fsa
for n in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do
    cat $INPUT >> $OUTDIR/temp.fsa
done

$VSEARCH --shuffle $OUTDIR/temp.fsa --output $OUTDIR/temp2.fsa

INPUT=$OUTDIR/temp2.fsa

echo 

echo Running vsearch

/usr/bin/time $VSEARCH \
    --derep_fulllength $INPUT \
    --threads $THREADS \
    --strand both \
    --minuniquesize $MINSIZE \
    --maxuniquesize $MAXSIZE \
    --topn $TOPN \
    --sizeout \
    --output $OUTDIR/derep.v.fsa \
    --uc $OUTDIR/derep.v.uc


echo

echo Running usearch

/usr/bin/time $USEARCH \
    --derep_fulllength $INPUT \
    --threads $THREADS \
    --strand both \
    --minuniquesize $MINSIZE \
    --topn $TOPN \
    --sizeout \
    --output $OUTDIR/derep.u.fsa \
    --uc $OUTDIR/derep.u.uc

