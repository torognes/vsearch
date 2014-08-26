#!/bin/sh

P=$1

DIR=.
DB=../data/Rfam_11_0.fasta
VSEARCH=../src/vsearch
USEARCH=$(which usearch)
SEED=1
THREADS=0
DUPLICATES=100

if [ "$P" = "u" ]; then
    PROG=$USEARCH
else
    if [ "$P" = "v" ]; then
	PROG=$VSEARCH
    else
	echo You must specify u or v as first argument
	exit
    fi
fi

echo Creating random test set

../src/vsearch --shuffle $DB --output $DIR/temp.fsa --seed $SEED > /dev/null 2> /dev/null
./select.pl $DIR/temp.fsa $DIR/q.fsa $DIR/db.fsa


cat q.fsa >> qq.fsa
for (( i=2; i <= $DUPLICATES; i++ )); do
    cat q.fsa >> qq.fsa
done

echo

echo Running search

/usr/bin/time $PROG \
    --usearch_global $DIR/qq.fsa \
    --db $DIR/db.fsa \
    --uc $DIR/uc.$P.txt \
    --blast6out $DIR/b6.$P.txt \
    --minseqlength 1 \
    --id 0.5 \
    --maxaccepts 1 \
    --maxrejects 32 \
    --strand plus \
    --threads $THREADS \
    --fulldp

echo

echo Results

./stats.pl $(grep -c "^>" $DIR/qq.fsa) $DIR/uc.$P.txt
