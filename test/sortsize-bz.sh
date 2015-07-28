#!/bin/bash

INPUT=../data/BioMarKs.fsa.bz2
TEMP=temp-bz.fsa

../src/vsearch-bz --shuffle $INPUT --output $TEMP

echo

/usr/bin/time ../src/vsearch-bz --sortbysize $TEMP --output vsortsize-bz.fsa \

echo

/usr/bin/time usearch --sortbysize $TEMP  -output usortsize.fsa
