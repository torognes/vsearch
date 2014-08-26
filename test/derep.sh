#!/bin/sh

INPUT=../data/PR2-18S-rRNA-V4.fsa

/usr/bin/time ../src/vsearch \
    --derep_fulllength $INPUT \
    --threads 1 \
    --strand both \
    --output vderep.fsa \
    --sizeout \
    --uc vderep.uc

/usr/bin/time usearch \
    --derep_fulllength $INPUT \
    --threads 1 \
    --strand both \
    --output uderep.fsa \
    --sizeout \
    --uc uderep.uc
