#!/bin/sh

INPUT=../data/BioMarKs.fsa
TEMP=temp.fsa

../src/vsearch --shuffle $INPUT --output $TEMP

echo

/usr/bin/time ../src/vsearch --sortbysize $TEMP --output vsortsize.fsa \

echo

/usr/bin/time usearch --sortbysize $TEMP  -output usortsize.fsa