#!/bin/sh

Q=../data/Rfam_9_1.fasta
DB=../data/Rfam_9_1.fasta
T=0
ID=0.5
MR=32

/usr/bin/time ../src/vsearch \
    --usearch_global $Q \
    --db $DB \
    --threads $T \
    --strand plus \
    --id $ID \
    --match 2 \
    --mismatch -4 \
    --gapopen 20I/2E \
    --gapext 2I/1E \
    --maxaccepts 1 \
    --maxrejects $MR \
    --alnout alnout.v.txt

/usr/bin/time usearch \
    --usearch_global $Q \
    --db $DB \
    --threads $T \
    --strand plus \
    --id $ID \
    --match 2 \
    --mismatch -4 \
    --gapopen 20I/2E \
    --gapext 2I/1E \
    --maxaccepts 1 \
    --maxrejects $MR \
    --alnout alnout.u.txt
