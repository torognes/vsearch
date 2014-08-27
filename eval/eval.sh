#!/bin/sh

P=$1

SEED=1
THREADS=0
DUPLICATES=10
DIR=.
DB=../data/Rfam_11_0.fasta

USEARCH=$(which usearch)

if [ $(uname -s) == "Linux" ]; then
    VSEARCH=$(ls -t ../bin/vsearch*linux* | head -1)
else
    VSEARCH=$(ls -t ../bin/vsearch*macosx* | head -1)
fi

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

echo Creating random test set

$VSEARCH --shuffle $DB --output $DIR/temp.fsa --seed $SEED > /dev/null 2> /dev/null
./select.pl $DIR/temp.fsa $DIR/q.fsa $DIR/db.fsa

cat q.fsa > qq.fsa
for (( i=2; i <= $DUPLICATES; i++ )); do
    cat q.fsa >> qq.fsa
done

echo

echo Running search

/usr/bin/time $PROG \
    --usearch_global $DIR/qq.fsa \
    --db $DIR/db.fsa \
    --id 0.5 \
    --maxaccepts 1 \
    --maxrejects 32 \
    --strand plus \
    --threads $THREADS \
    --userout $DIR/userout.$P.txt \
    --userfields query+target+id+qcov

#    --fulldp
#    --minseqlength 1

echo

echo Results

./stats.pl $(grep -c "^>" $DIR/qq.fsa) $DIR/userout.$P.txt

rm temp.fsa q.fsa qq.fsa db.fsa userout.$P.txt
