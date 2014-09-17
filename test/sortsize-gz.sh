#!/bin/sh

INPUT=../data/BioMarKs.fsa.gz
TEMP=temp-gz.fsa

../src/vsearch-gz --shuffle $INPUT --output $TEMP

echo

/usr/bin/time ../src/vsearch-gz --sortbysize $TEMP --output vsortsize-gz.fsa \

echo

/usr/bin/time usearch --sortbysize $TEMP  -output usortsize-gz.fsa
