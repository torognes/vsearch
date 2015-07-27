#!/bin/sh

INPUT=../data/BioMarKs.fsa

echo

/usr/bin/time ../src/vsearch --fastx_subsample $INPUT --fastaout subsampled.fsa --sample_pct 10 --sizein --sizeout

