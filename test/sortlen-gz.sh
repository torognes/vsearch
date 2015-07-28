#!/bin/bash

INPUT=../data/PR2-18S-rRNA-V4.fsa.gz

echo

/usr/bin/time ../src/vsearch-gz --sortbylength $INPUT --output vsortlen-gz.fsa

echo

/usr/bin/time usearch --sortbylength $INPUT --output usortlen-gz.fsa
