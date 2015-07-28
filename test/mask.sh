#!/bin/bash

DB=../data/BioMarKs50k.fsa
VSEARCH=../src/vsearch
USEARCH=$(which usearch)

MASK=dust
WIDTH=80

/usr/bin/time $VSEARCH \
    --maskfasta $DB \
    --fasta_width $WIDTH \
    --qmask $MASK \
    --hardmask \
    --output masked.v.fsa

/usr/bin/time $USEARCH \
    --maskfasta $DB \
    --qmask $MASK \
    --hardmask \
    --output masked.u.fsa

