#!/bin/sh

DB=../data/BioMarKs.fsa
VSEARCH=../src/vsearch
USEARCH=$(which usearch)

/usr/bin/time $VSEARCH \
    --maskfasta $DB \
    --qmask dust \
    --hardmask \
    --output masked.v.fsa

/usr/bin/time $USEARCH \
    --maskfasta $DB \
    --qmask dust \
    --hardmask \
    --output masked.u.fsa
