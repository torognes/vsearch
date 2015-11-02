#!/bin/bash

DATA=../../vsearch-data/BioMarKs.fsa
INPUT=temp1.fsa
OUTPUT=temp2.fsa
VSEARCH=../bin/vsearch

echo Generating a dataset

rm -f $INPUT

cat $DATA >> $INPUT
cat $DATA >> $INPUT
cat $DATA >> $INPUT
cat $DATA >> $INPUT
cat $DATA >> $INPUT
cat $DATA >> $INPUT
cat $DATA >> $INPUT
cat $DATA >> $INPUT
cat $DATA >> $INPUT
cat $DATA >> $INPUT

echo Shuffling

/usr/bin/time $VSEARCH --shuffle $INPUT --output $OUTPUT

rm $INPUT $OUTPUT
