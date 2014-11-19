#!/bin/bash

DB=../data/simm/simm_sp.fa

THREADS=0

UCHIME=$(which uchime)
USEARCH=$(which usearch)
VSEARCH=../src/vsearch

echo Evaluation of uchime_ref

for t in - m1 m2 m3 m4 m5 i1 i2 i3 i4 i5; do

    echo
    echo Running programs on dataset $t
    echo

    if [ "$t" == "-" ]; then
        INPUT=../data/simm/simm.fa
    else
        INPUT=../data/simm/simm.$t.fa
    fi

    $UCHIME --input $INPUT --db $DB --uchimeout o.$t.uchimeout --minh 0.28 --mindiv 0.8 ; grep Y$ o.$t.uchimeout | cut -f2 > o.$t.chimeras

    $VSEARCH --uchime_ref $INPUT --db $DB --chimeras v.$t.chimeras

    $USEARCH --uchime_ref $INPUT --db $DB --strand plus --chimeras u.$t.chimeras

done

echo
echo -e  "\t__________m=2__________\t__________m=3__________\t__________m=4__________"
echo -ne "Div/Evo\tUSEARCH\tUCHIME\tVSEARCH"
echo -ne "\tUSEARCH\tUCHIME\tVSEARCH"
echo -e  "\tUSEARCH\tUCHIME\tVSEARCH"

for r in 97_99 95_97 90_95; do
    for t in - i1 i2 i3 i4 i5 m1 m2 m3 m4 m5; do
        EXT="$t.chimeras"
        echo -ne "$r$t"
        for m in 2 3 4; do
            echo -ne "\t$(grep -c _m${m}_${r} u.$EXT)"
            echo -ne "\t$(grep -c _m${m}_${r} o.$EXT)"
            echo -ne "\t$(grep -c _m${m}_${r} v.$EXT)"
        done
        echo
    done
    echo
done

rm -f o.* u.* v.*

echo Done
