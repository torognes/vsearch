#!/bin/bash

P=$1

# Searching against dereplicated sequences with 
# very fasta names
#Q=../../vsearch-data/PR2-18S-rRNA-V4.fsa
#DB=../../vsearch-data/PR2-18S-rRNA-V4.derep.fsa

# Searching with a compressed Q and when
# Q and DB are not identical
#Q=../../vsearch-data/BioMarKs.fsa.gz
#DB=../../vsearch-data/BioMarKs50k.fsa

# Search with longer reads
Q=../../vsearch-data/constaint.fsa.bz2
DB=../../vsearch-data/constaint.fsa.bz2


USEARCH=$(which usearch)
VSEARCH=../bin/vsearch

if [ "$P" == "u" ]; then
    PROG=$USEARCH
else
    if [ "$P" == "v" ]; then
        PROG=$VSEARCH
    else
        echo You must specify u or v as first argument
        exit
    fi
fi

CMD="/usr/bin/time $PROG \
    --search_exact $Q \
    --db $DB \
    --strand plus \
    --sizeout \
    --alnout alnout.$P.txt \
    --samout samout.$P.txt \
    --uc uc.$P.fsa \
    --blast6out blast6out.$P.bl6 \
    --dbmatched dbmatched.$P.fsa \
    --dbnotmatched dbnotmatched.$P.fsa \
    --matched matched.$P.fsa \
    --notmatched notmatched.$P.fsa"

echo Search test
echo
echo Running command: $CMD
echo

$CMD
