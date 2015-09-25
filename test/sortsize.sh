#!/bin/bash

INPUT=../../vsearch-data/BioMarKs.fsa
TEMP=temp.fsa

../bin/vsearch --shuffle $INPUT --output $TEMP

echo

/usr/bin/time ../bin/vsearch --sortbysize $TEMP --output vsortsize.fsa \

echo

/usr/bin/time usearch --sortbysize $TEMP  -output usortsize.fsa