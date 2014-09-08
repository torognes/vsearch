#!/bin/sh

INPUT=../data/PR2-18S-rRNA-V4.fsa.bz2

/usr/bin/time ../src/vsearch-bz \
    --derep_fulllength $INPUT \
    --threads 1 \
    --strand both \
    --output vderep-bz.fsa \
    --sizeout \
    --uc vderep-bz.uc

/usr/bin/time usearch \
    --derep_fulllength $INPUT \
    --threads 1 \
    --strand both \
    --output uderep-bz.fsa \
    --sizeout \
    --uc uderep-bz.uc
