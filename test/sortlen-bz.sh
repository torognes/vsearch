#!/bin/sh

INPUT=../data/PR2-18S-rRNA-V4.fsa.bz2

echo

/usr/bin/time ../src/vsearch-bz --sortbylength $INPUT --output vsortlen-bz.fsa

echo

/usr/bin/time usearch --sortbylength $INPUT --output usortlen.fsa
