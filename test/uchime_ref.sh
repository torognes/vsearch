#!/bin/sh

cd ../test

P=$1

INPUT=../data/simm/simm.m1.fa
DB=../data/simm/simm_sp.fa

#INPUT=../data/PR2-18S-rRNA-V4.derep.fsa
#DB=../data/PR2-18S-rRNA-V4.ref.fsa

THREADS=0

UCHIME=$(which uchime)
USEARCH=$(which usearch)
VSEARCH=../src/vsearch

if [ "$P" == "u" ]; then
    PROG=$USEARCH
else
    if [ "$P" == "v" ]; then
        PROG=$VSEARCH
    else
        if [ "$P" == "o" ]; then
            PROG=$UCHIME
        else
            echo You must specify u or v or o as first argument
            exit
        fi
    fi
fi

if [ "$P" == "o" ]; then

    CMD="/usr/bin/time $PROG \
        --input $INPUT \
        --db $DB \
        --uchimeout $P.uchimeout \
        --uchimealns $P.uchimealns \
        --minh 0.28 \
        --mindiv 0.8"

else

    CMD="/usr/bin/time $PROG \
      --uchime_ref $INPUT \
      --db $DB \
      --strand plus \
      --chimeras $P.chimeras \
      --nonchimeras $P.nonchimeras \
      --uchimealns $P.uchimealns \
      --uchimeout $P.uchimeout"

fi
    
echo uchime_ref test
echo
echo Running command: $CMD
echo

$CMD
