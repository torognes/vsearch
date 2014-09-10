#!/bin/sh

INPUT=../data/PR2-18S-rRNA-V4.fsa.gz

/usr/bin/time ../src/vsearch-gz \
    --derep_fulllength $INPUT \
    --threads 1 \
    --strand both \
    --output vderep-gz.fsa \
    --sizeout \
    --uc vderep-gz.uc

/usr/bin/time usearch \
    --derep_fulllength $INPUT \
    --threads 1 \
    --strand both \
    --output uderep-gz.fsa \
    --sizeout \
    --uc uderep-gz.uc
