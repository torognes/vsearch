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


exit 0
