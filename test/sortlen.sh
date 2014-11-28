#!/bin/sh

INPUT=../data/PR2-18S-rRNA-V4.fsa

echo

/usr/bin/time ../src/vsearch --sortbylength $INPUT --output vsortlen.fsa --sizeout

echo

/usr/bin/time usearch --sortbylength $INPUT --output usortlen.fsa --sizeout
