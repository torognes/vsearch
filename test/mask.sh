#!/bin/sh

DB=../data/BioMarKs.fsa

/usr/bin/time ../src/vsearch \
    --mask $DB \
    --dbmask dust \
    --hardmask \
    --minseqlen 1 \
    --fasta_width 60 \
    --output masked.fsa
