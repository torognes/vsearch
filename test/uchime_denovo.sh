#!/bin/bash

cd ../test

P=$1

INPUT=../data/PR2-18S-rRNA-V4.derep.fsa

THREADS=0

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
    --uchime_denovo $INPUT \
    --chimeras $P.chimeras \
    --nonchimeras $P.nonchimeras \
    --uchimealns $P.uchimealns \
    --uchimeout $P.uchimeout"

echo uchime_denovo test
echo
echo Running command: $CMD
echo

$CMD
