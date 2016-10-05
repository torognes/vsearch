#!/bin/bash -

## Print a header
SCRIPT_NAME="Unclassified tests"
LINE=$(printf "%076s\n" | tr " " "-")
printf "# %s %s\n" "${LINE:${#SCRIPT_NAME}}" "${SCRIPT_NAME}"

## Declare a color code for test results
RED="\033[1;31m"
GREEN="\033[1;32m"
NO_COLOR="\033[0m"

failure () {
    printf "${RED}FAIL${NO_COLOR}: ${1}\n"
    # exit -1
}

success () {
    printf "${GREEN}PASS${NO_COLOR}: ${1}\n"
}


## Is vsearch installed?
VSEARCH=$(which vsearch)
DESCRIPTION="check if vsearch is in the PATH"
[[ "${VSEARCH}" ]] && success "${DESCRIPTION}" || failure "${DESCRIPTION}"

#*****************************************************************************#
#                                                                             #
#                    Clustering UC format CIGAR alignment                     #
#                                                                             #
#*****************************************************************************#

## usearch 6, 7 and 8 output a "=" when the sequences are identical
DESCRIPTION="CIGAR alignment is \"=\" when the sequences are identical"
UC_OUT=$("${VSEARCH}" \
             --cluster_fast <(printf ">seq1\nACGT\n>seq2\nACGT\n") \
             --id 0.97 \
             --quiet \
             --minseqlength 1 \
             --uc - | grep "^H" | cut -f 8)

[[ "${UC_OUT}" == "=" ]] && \
    success  "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

## clean
unset UC_OUT


## is the 3rd column of H the query length or the alignment length?
DESCRIPTION="3rd column of H is the query length"
UC_OUT=$("${VSEARCH}" \
             --cluster_fast <(printf ">seq1\nACGT\n>seq2\nACAGT\n") \
             --id 0.5 \
             --quiet \
             --minseqlength 1 \
             --uc - | grep "^H")

awk 'BEGIN {FS = "\t"} {$3 == 4 && $9 == "seq1"}' <<< "${UC_OUT}" && \
    success  "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

## clean
unset UC_OUT


#*****************************************************************************#
#                                                                             #
#                        UC format when dereplicating                         #
#                                                                             #
#*****************************************************************************#

## sizein is taken into account
DESCRIPTION="when prefix dereplicating, --uc output accounts for --sizein"
s=$(printf ">seq1;size=3;\nACGT\n>seq2;size=1;\nACGT\n" | \
           "${VSEARCH}" \
               --derep_prefix - \
               --quiet \
               --sizein \
               --minseqlength 1 \
               --uc - | grep "^C" | cut -f 3)

(( ${s} == 4 )) && \
    success  "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

# clean
unset s


## vsearch reports H record when sequences have the same length
DESCRIPTION="when prefix dereplicating same length sequences, --uc reports H record"
H=$(printf ">seq1\nACGT\n>seq2\nACGT\n" | \
           "${VSEARCH}" \
               --derep_prefix - \
               --quiet \
               --minseqlength 1 \
               --uc - | grep "^H")

[[ -n ${H} ]] && \
    success  "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

# clean
unset H

## vsearch reports H record when sequences have different lengths
DESCRIPTION="when prefix dereplicating a shorter sequence, --uc reports H record"
H=$(printf ">seq1\nACGTA\n>seq2\nACGT\n" | \
           "${VSEARCH}" \
               --derep_prefix - \
               --quiet \
               --minseqlength 1 \
               --uc - | grep "^H")

[[ -n ${H} ]] && \
    success  "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

## clean
unset H


## --derep_prefix does not support the option --strand
DESCRIPTION="--derep_prefix does not support the option --strand"
printf ">seq1\nAATT\n>seq2\nTTAA\n" | \
    "${VSEARCH}" \
        --derep_prefix - \
        --quiet \
        --strand both \
        --minseqlength 1 \
        --uc - 2> /dev/null && \
    failure "${DESCRIPTION}" || \
        success  "${DESCRIPTION}"

# clean
unset H


## --derep_fulllength accepts the option --strand
DESCRIPTION="--derep_fulllength accepts the option --strand"
printf ">seq1\nAATT\n>seq2\nTTAA\n" | \
    "${VSEARCH}" \
        --derep_fulllength - \
        --quiet \
        --strand both \
        --minseqlength 1 \
        --uc /dev/null 2> /dev/null && \
    success  "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

# clean
unset H


## --derep_fulllength searches both strands
DESCRIPTION="--derep_fulllength searches both strands"
C=$(printf ">seq1\nAACC\n>seq2\nGGTT\n" | \
           "${VSEARCH}" \
               --derep_fulllength - \
               --quiet \
               --strand both \
               --minseqlength 1 \
               --uc - | grep -c "^C")

# There should be only cluster
(( ${C} == 1 )) && \
    success  "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

# clean
unset C

#*****************************************************************************#
#                                                                             #
#         fastq_trunclen and discarded short sequences (issue 203)            #
#                                                                             #
#*****************************************************************************#

DESCRIPTION="entries shorter than the --fastq_trunclength value are discarded"
"${VSEARCH}" \
    --fastq_filter <(printf "@seq1\nACGT\n+\nIIII\n") \
    --fastq_trunclen 5 \
    --quiet \
    --fastqout - \
    2> /dev/null | \
    grep -q "seq1" && \
    failure "${DESCRIPTION}" || \
        success  "${DESCRIPTION}"

DESCRIPTION="entries equal or longer than the --fastq_trunclength value are kept"
"${VSEARCH}" \
    --fastq_filter <(printf "@seq1\nACGT\n+\nIIII\n") \
    --fastq_trunclen 4 \
    --quiet \
    --fastqout - \
    2> /dev/null | \
    grep -q "seq1" && \
    success  "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"


#*****************************************************************************#
#                                                                             #
#     fastx_filter ignores sizein when relabeling fasta input (issue #204)    #
#                                                                             #
#*****************************************************************************#

# https://github.com/torognes/vsearch/issues/204
#
# --fastx_filter ignores input sequence abundances when relabeling
# with fasta input, --sizein and --sizeout options
DESCRIPTION="fastx_filter reports sizein when relabeling fasta (issue #204)"
"${VSEARCH}" \
    --fastx_filter <(printf ">seq1;size=5;\nACGT\n") \
    --sizein \
    --relabel_md5 \
    --sizeout \
    --quiet \
    --fastaout - \
    2> /dev/null | \
    grep -q ";size=5;" && \
    success  "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"


#*****************************************************************************#
#                                                                             #
#         check pairwise alignment correctness (Flouri et al., 2015)          #
#                                                                             #
#*****************************************************************************#

# http://biorxiv.org/content/early/2015/11/12/031500

# In USEARCH and VSEARCH the gap opening penalty includes the gap
# extension penalty of the first residue, while in other programs it
# does not. So if the gap open penalty is 40 and the gap extension
# penalty is 1, then a single nucleotide gap will get a penalty of 40
# in USEARCH and VSEARCH, and 41 in other programs.

# In Flouri's tests, the gap opening penalty does not include the gap
# extension penalty, and the optimal alignments contain two
# independent gaps. Therefore, USEARCH and VSEARCH should return score
# values equal to the scores indicated by Flouri, minus twice the gap
# extension penalty (e.g., a score of -72 reported by Flouri
# corresponds to a score of -70 with USEARCH and VSEARCH). The
# expected score values in the tests below take that into account.

# test 1 requires the possibility to set independent match/mismatch
# scores for the different pairs of nucleotides. Not possible to
# replicate in vsearch: ">seq1\nGGTGTGA\n>seq2\nTCGCGT\n"

# test 2 uses a match score of zero, not possible with vsearch (Fatal
# error: The argument to --match must be positive)
# ">seq1\nAAAGGG\n>seq2\nTTAAAAGGGGTT\n"

# test 3 (score should be -70 in USEARCH/VSEARCH)
DESCRIPTION="Flouri 2015 pairwise alignment correctness tests (test 3)"
score=$("${VSEARCH}" \
            --allpairs_global <(printf ">seq1\nAAATTTGC\n>seq2\nCGCCTTAC\n") \
            --acceptall \
            --gapopen 40 \
            --gapext 1\
            --match 10 \
            --mismatch -30 \
            --qmask none \
            --quiet \
            --userfields raw \
            --userout -)

(( ${score} == -70 )) && \
    success  "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

# test 4 (score should be -60 in USEARCH/VSEARCH
DESCRIPTION="Flouri 2015 pairwise alignment correctness tests (test 4)"
score=$("${VSEARCH}" \
            --allpairs_global <(printf ">seq1\nTAAATTTGC\n>seq2\nTCGCCTTAC\n") \
            --acceptall \
            --gapopen 40 \
            --gapext 1\
            --match 10 \
            --mismatch -30 \
            --qmask none \
            --quiet \
            --userfields raw \
            --userout -)

(( ${score} == -60 )) && \
    success  "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

# test 5 (identical to test 3)

# test 6 (score should be -44 in USEARCH/VSEARCH
DESCRIPTION="Flouri 2015 pairwise alignment correctness tests (test 6)"
score=$("${VSEARCH}" \
            --allpairs_global <(printf ">seq1\nAGAT\n>seq2\nCTCT\n") \
            --acceptall \
            --gapopen 25 \
            --gapext 1\
            --match 10 \
            --mismatch -30 \
            --qmask none \
            --quiet \
            --userfields raw \
            --userout -)

(( ${score} == -44 )) && \
    success  "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

exit 0
