#!/bin/bash -
# Byte-identity matrix between two vsearch binaries.
#
# Verification step 2 of TBD_20260804_c_style_elimination.md: replacing a
# format string with a typed writer must not change a single output byte, so
# every commit of that migration is checked by running the same command
# matrix against a baseline binary and this branch's binary and comparing
# every file both produce.
#
# usage: bash byte_identity.sh <binary_a> <binary_b> [--threads N]
#
# Ported from swarm's byte_identity.sh, with vsearch's per-command option
# matrix. The default is --threads 1: with more than one thread the record
# order of the search commands is not deterministic, and every difference
# reported would be reordering rather than a formatting change. Pass
# --threads 8 to exercise the interleaving hazard instead (a split output
# line is several stdio calls where it used to be one); that run is compared
# after sorting.
#
# The --log output is compared with its volatile lines removed (the command
# line, wall-clock timings, memory figures and the host name).

BIN_A="${1}"
BIN_B="${2}"
THREADS=1
[[ "${3}" == "--threads" ]] && THREADS="${4}"

[[ -x "${BIN_A}" ]] || { echo "not executable: ${BIN_A}" ; exit 1 ; }
[[ -x "${BIN_B}" ]] || { echo "not executable: ${BIN_B}" ; exit 1 ; }

# each run happens inside its own output directory, so a relative binary
# path would no longer resolve there
BIN_A=$(cd "$(dirname "${BIN_A}")" && pwd)/$(basename "${BIN_A}")
BIN_B=$(cd "$(dirname "${BIN_B}")" && pwd)/$(basename "${BIN_B}")

# Compare like with like. A DEBUG build carries -D_GLIBCXX_DEBUG and full
# debug info and is several times the size of a RELEASE build; the two also
# differ in what their memory figures and their assertion behaviour report,
# so a DEBUG-vs-RELEASE run reports differences that have nothing to do with
# the code under test -- and reports them intermittently, which is worse
# than reporting them always.
SIZE_A=$(stat -c %s "${BIN_A}")
SIZE_B=$(stat -c %s "${BIN_B}")
if [[ $(( SIZE_A > SIZE_B ? SIZE_A / SIZE_B : SIZE_B / SIZE_A )) -ge 2 ]] ; then
    echo "REFUSING: ${SIZE_A} vs ${SIZE_B} bytes -- these look like different build types"
    echo "  (${BIN_A})"
    echo "  (${BIN_B})"
    echo "Rebuild both the same way (both --enable-debug, or both release) and re-run."
    exit 2
fi

WORKDIR=$(mktemp -d)
trap 'rm -rf "${WORKDIR}"' EXIT

# --- inputs -----------------------------------------------------------------
#
# Deliberately varied, because a formatting change shows up in the fields a
# plain input never populates: abundances that tie (the label tie-break),
# a header carrying a byte >= 0x80 (element_order), a header with a space
# (label truncation), sequences that are prefixes of one another
# (--derep_prefix), ambiguities and lowercase (--qmask, --fastx_mask),
# lengths that are and are not a multiple of the FASTA wrap width, and a
# sequence pair far enough apart to produce a non-trivial CIGAR.

QUERY="${WORKDIR}/query.fas"
cat > "${QUERY}" << 'END_OF_FASTA'
>q1;size=5;
ACGTACGTACGTAAAACCCGGGTTTACGTACGTACGTAAAACCCGGGTTTACGTACGTACGTAAAACCCGGGTTTAAAA
>q2;size=5;
ACGTACGTACGTAAAACCCGGGTTTACGTACGTACGTAAAACCCGGGTTTACGTACGTACGTAAAACCCGGGTTTAAA
>q3 with a space;size=3;
ACGTACGTACGTAAAACCCGGGTTTACGTACGTACGTAAAACCCGGGTTTACGTACGTACGTAAAACCCGGGTTTAAAAC
>q4sø;size=1;
ACGTACGTACGTAAAACCCGGGTTTACGTACGTACGTAAAACCCGGGTTTACGTACGTACGTAAAACCCGGGTTTAAAT
>q5;size=1;
acgtacgtacgtaaaacccgggtttacgtacgtacgtaaaacccgggtttacgtacgtacgtaaaacccgggtttaata
>q6;size=1;
ACGTACGTACGTAAAACCCGGGTTTACGTNCGTACGTAAAACCCGGGTTTACGTACGTACGTAAAACCCGGGTTTAATA
>q7;size=2;
TTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAA
>q8;size=1;
TTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAG
END_OF_FASTA

# reference database: q1 and q7 exactly, plus one sequence no query matches
# and one that only matches on the minus strand
DB="${WORKDIR}/db.fas"
cat > "${DB}" << 'END_OF_FASTA'
>ref1;size=10;
ACGTACGTACGTAAAACCCGGGTTTACGTACGTACGTAAAACCCGGGTTTACGTACGTACGTAAAACCCGGGTTTAAAA
>ref2;size=4;
TTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAA
>ref3;size=1;
GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
>ref4revcomp;size=1;
TTTTAAACCCGGGTTTTGCTACCGTACGTACGTAAACCCGGGTTTTGCTACCGTACGTACGTAAACCCGGGTTTACGT
END_OF_FASTA

# taxonomy-annotated reference, for --sintax
TAXDB="${WORKDIR}/taxdb.fas"
cat > "${TAXDB}" << 'END_OF_FASTA'
>ref1;tax=d:Bacteria,p:Firmicutes,c:Bacilli,o:Bacillales,f:Bacillaceae,g:Bacillus;
ACGTACGTACGTAAAACCCGGGTTTACGTACGTACGTAAAACCCGGGTTTACGTACGTACGTAAAACCCGGGTTTAAAA
>ref2;tax=d:Bacteria,p:Proteobacteria,c:Gammaproteobacteria,o:Enterobacterales,f:Enterobacteriaceae,g:Escherichia;
TTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAA
>ref3;tax=d:Archaea,p:Euryarchaeota;
GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
END_OF_FASTA

# prefix-nested input, for --derep_prefix
PREFIX="${WORKDIR}/prefix.fas"
cat > "${PREFIX}" << 'END_OF_FASTA'
>p1;size=3;
ACGTACGTACGTAAAACCCGGGTTTACGTACGTACGT
>p2;size=2;
ACGTACGTACGTAAAACCCGGGTTT
>p3;size=1;
ACGTACGTACGTAAAA
>p4;size=1;
TTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTG
END_OF_FASTA

# a chimera: the first half of ref1 followed by the second half of ref2,
# at a low abundance so that --abskew accepts it
CHIM="${WORKDIR}/chim.fas"
cat > "${CHIM}" << 'END_OF_FASTA'
>parentA;size=50;
ACGTACGTACGTAAAACCCGGGTTTACGTACGTACGTAAAACCCGGGTTTACGTACGTACGTAAAACCCGGGTTTAAAAACGTACGTACGTAAAACCCGGGTTTACGTACGTACGTAAAACCCGGGTTT
>parentB;size=40;
TTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAA
>achimera;size=1;
ACGTACGTACGTAAAACCCGGGTTTACGTACGTACGTAAAACCCGGGTTTACGTACGTACGTAAAACCCGGGTTTAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAA
>another;size=1;
ACGTACGTACGTAAAACCCGGGTTTACGTACGTACGTAAAACCCGGGTTTACGTACGTACGTAAAACCCGGGTTTAAAAACGTACGTACGTAAAACCCGGGTTTACGTACGTACGTAAAACCCGGGTTA
END_OF_FASTA

# FASTQ inputs. The quality line is generated from the sequence length rather
# than written out, so the two cannot drift apart; the phred values cycle over
# a wide range so that every --fastq_* statistic (min, max, mean, expected
# errors, the eestats precision ladder) has something to report.
make_fastq() {
    # <fasta in> <fastq out> [starting phred offset]
    awk -v seed="${3:-0}" '
        /^>/ { header = substr($0, 2); next }
        {
          quality = ""
          for (position = 1; position <= length($0); ++position) {
            quality = quality sprintf("%c", 33 + ((position * 7 + seed) % 41))
          }
          print "@" header; print $0; print "+"; print quality
          seed += 13
        }' "${1}" > "${2}"
}

FASTQ="${WORKDIR}/reads.fastq"
make_fastq "${QUERY}" "${FASTQ}"

# A paired set for --fastq_mergepairs / --fastq_join / --fastx_syncpairs. R2 is
# the reverse complement of the tail of the same fragment, so the pair overlaps
# by 2 * read length - fragment length; one pair carries a mismatch inside the
# overlap (so the merged quality recomputation runs) and one fragment is long
# enough that its reads do not overlap at all (so the not-merged outputs are
# populated too).
revcomp() { printf '%s' "${1}" | rev | tr 'ACGTacgt' 'TGCAtgca' ; }

FRAG_A="ACGTACGTACGTAAAACCCGGGTTTACGTACGTACGTAAAACCCGGGTTTACGTACGTACGTAAAACCCGGGTTTAAAACCCGGGTTTACGTACGTACGT"
FRAG_B="TTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTT"
FRAG_C="ACGTACGTACGTAAAACCCGGGTTTACGTACGTACGTAAAACCCGGGTTTACGTACGTACGTAAAACCCGGGTTTAAAACCCGGGTTTACGTACGTACGTAAAACCCGGGTTTACGTACGTACGTAAAACCCGGGTTTACGTACGTACGT"

FASTA_R1="${WORKDIR}/r1.fas"
FASTA_R2="${WORKDIR}/r2.fas"
{
    printf '>p1\n%s\n' "${FRAG_A:0:70}"
    printf '>p2\n%s\n' "${FRAG_B:0:70}"
    printf '>p3\n%s\n' "${FRAG_C:0:70}"
} > "${FASTA_R1}"
{
    printf '>p1\n%s\n' "$(revcomp "${FRAG_A: -70}")"
    # one mismatch inside the overlap: the leading base of the R2 read differs
    printf '>p2\n%s\n' "G$(revcomp "${FRAG_B: -70}" | cut -c2-)"
    printf '>p3\n%s\n' "$(revcomp "${FRAG_C: -70}")"
} > "${FASTA_R2}"

R1="${WORKDIR}/r1.fastq"
R2="${WORKDIR}/r2.fastq"
make_fastq "${FASTA_R1}" "${R1}" 0
make_fastq "${FASTA_R2}" "${R2}" 5

LABELS="${WORKDIR}/labels.txt"
printf 'ref1\nref3\n' > "${LABELS}"

# The UDB is built once, by BIN_A, and read by both: a byte-identity run
# compares readers and writers, not two independently built databases.
UDB="${WORKDIR}/db.udb"
"${BIN_A}" --makeudb_usearch "${DB}" --output "${UDB}" --quiet 2> /dev/null

# --- the run harness --------------------------------------------------------

# Output files a run may produce. Every one is compared; the ones a given
# command does not write are absent in both directories, which counts as
# agreement (see compare_run).
COMPARED=(stdout stderr status log out.fas out.fastq out.txt out.uc out.b6
          out.aln out.user out.sam out.pairs out.tab out.otu out.biom
          out.mothur out.msa out.cons out.prof out.cent out.chim out.nonchim
          out.borderline out.alns out.uchimeout out.notmatched out.rejected
          out.r1 out.r2 out.merged out.notmerged1 out.notmerged2 out.eestats)

differences=0
comparisons=0
runs=0
failed_runs=0
FAILED_LABELS=()

# Two things in an output file legitimately differ between the two binaries
# and are not formatting: the path each one was invoked as (echoed in the log,
# in --alnout's first line and in --help's usage line) and, in the log, the
# wall-clock timings and memory figures. Everything else -- and the log is
# itself a converted writer -- is stable between two builds of the same tree.
normalize() {
    sed -E -e "s@${BIN_A}@BINARY@g" -e "s@${BIN_B}@BINARY@g" \
           -e 's/[0-9]+\.[0-9]+ ?(seconds|sec)\b/TIME/g' \
           -e '/^Started  /d' \
           -e '/^Finished /d' \
           -e '/[Ee]lapsed time/d' \
           -e '/[Mm]aximum resident/d' \
           -e '/[Rr]ate:/d' \
           -e '/^Max memory /d' \
           "${1}"
}

# --cut is the one command in the matrix with no threaded phase, and vsearch
# rejects an option a command does not accept rather than ignoring it, so
# --threads cannot be passed unconditionally.
rejects_threads() { [[ "${1}" == "--cut" ]] ; }

run_one() {
    # $1 = binary, $2 = output directory, rest = options
    #
    # Run inside the output directory with relative output names: vsearch
    # echoes each output path in its log, so absolute paths carrying the
    # directory name would make the log differ for that reason alone.
    local binary="${1}" outdir="${2}"
    shift 2
    local threads=(--threads "${THREADS}")
    rejects_threads "${1}" && threads=()
    mkdir -p "${outdir}"
    (
        cd "${outdir}" || exit 1
        "${binary}" "${@}" "${threads[@]}" --log log \
            > stdout 2> stderr
        echo "$?" > status
    )
}

compare_run() {
    local label="${1}"
    shift
    rm -rf "${WORKDIR}/a" "${WORKDIR}/b"
    run_one "${BIN_A}" "${WORKDIR}/a" "${@}"
    run_one "${BIN_B}" "${WORKDIR}/b" "${@}"
    runs=$(( runs + 2 ))
    local status
    for status in "$(cat "${WORKDIR}/a/status")" "$(cat "${WORKDIR}/b/status")" ; do
        if [[ "${status}" -ne 0 ]] ; then
            failed_runs=$(( failed_runs + 1 ))
            FAILED_LABELS+=("${label}")
        fi
    done
    local suffix file_a file_b
    for suffix in "${COMPARED[@]}" ; do
        file_a="${WORKDIR}/a/${suffix}"
        file_b="${WORKDIR}/b/${suffix}"
        # A command that does not write a given output writes it in neither
        # directory: both absent is agreement, one absent is not.
        [[ ! -e "${file_a}" && ! -e "${file_b}" ]] && continue
        comparisons=$(( comparisons + 1 ))
        normalize "${file_a}" > "${file_a}.clean"
        normalize "${file_b}" > "${file_b}.clean"
        file_a="${file_a}.clean"
        file_b="${file_b}.clean"
        if [[ "${THREADS}" -ne 1 ]] ; then
            # More than one thread means the record order is not fixed; the
            # interleaving hazard this mode exists for corrupts lines, which
            # a sorted comparison still catches.
            sort "${file_a}" -o "${file_a}" 2> /dev/null
            sort "${file_b}" -o "${file_b}" 2> /dev/null
        fi
        if ! cmp -s "${file_a}" "${file_b}" ; then
            differences=$(( differences + 1 ))
            echo "DIFF: ${label} [${suffix}]"
            diff "${file_a}" "${file_b}" | head -8
        fi
    done
}

# An input vsearch cannot read makes this whole script lie: both binaries
# abort in the parser, write the same error and the same status, and the
# matrix scores hundreds of identical comparisons without having compared
# any output at all. Check first, and say so loudly.
smoke_test() {
    local dir="${WORKDIR}/smoke"
    rm -rf "${dir}"
    run_one "${BIN_A}" "${dir}" --usearch_global "${QUERY}" --db "${DB}" \
            --id 0.9 --blast6out out.b6 --quiet
    local status
    status=$(cat "${dir}/status")
    if [[ "${status}" -ne 0 ]] || [[ ! -s "${dir}/out.b6" ]] ; then
        echo "REFUSING: no hits came out of the smoke test"
        echo "  status: ${status}"
        [[ -s "${dir}/stderr" ]] && echo "  stderr: $(tail -2 "${dir}/stderr")"
        echo "Every run below would fail the same way in both binaries, and the"
        echo "matrix would report BYTE-IDENTICAL without comparing anything."
        exit 2
    fi
    if [[ ! -s "${UDB}" ]] ; then
        echo "REFUSING: the UDB was not built; --udbinfo/--udbstats would compare nothing"
        exit 2
    fi
    rm -rf "${dir}"
}
smoke_test

# --- the matrix -------------------------------------------------------------
#
# One entry per output writer the migration touches, grouped by the phase of
# TBD_20260804_c_style_elimination.md that converts it.

# search: results.cpp (every output format), showalign, usearch_global,
# search_exact, dbindex, progress
for id_value in 0.75 0.97 1.0 ; do
    compare_run "usearch_global --id ${id_value}" \
        --usearch_global "${QUERY}" --db "${DB}" --id "${id_value}" \
        --alnout out.aln --blast6out out.b6 --uc out.uc --samout out.sam \
        --fastapairs out.pairs --matched out.fas --notmatched out.notmatched \
        --userout out.user --userfields 'query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+evalue+bits+ql+tl+qs+ts+pv+qcov+tcov+ids+gaps+qrow+trow+qstrand+tstrand+aln+caln+qilo+qihi+tilo+tihi+id0+id1+id2+id3+id4+exts+raw+qframe+tframe'
done
compare_run "usearch_global --strand both" \
    --usearch_global "${QUERY}" --db "${DB}" --id 0.8 --strand both \
    --alnout out.aln --uc out.uc --blast6out out.b6
compare_run "usearch_global --output_no_hits" \
    --usearch_global "${QUERY}" --db "${DB}" --id 0.99 --output_no_hits \
    --blast6out out.b6 --uc out.uc --userout out.user --userfields 'query+target+id'
compare_run "usearch_global --maxaccepts 0 --maxrejects 0" \
    --usearch_global "${QUERY}" --db "${DB}" --id 0.7 --maxaccepts 0 \
    --maxrejects 0 --alnout out.aln --uc out.uc
compare_run "usearch_global --db udb" \
    --usearch_global "${QUERY}" --db "${UDB}" --id 0.8 --alnout out.aln --uc out.uc
# --rowlen drives showalign.cpp's run-time field widths (phase 6)
compare_run "usearch_global --rowlen 0" \
    --usearch_global "${QUERY}" --db "${DB}" --id 0.8 --alnout out.aln --rowlen 0
compare_run "usearch_global --rowlen 17" \
    --usearch_global "${QUERY}" --db "${DB}" --id 0.8 --alnout out.aln --rowlen 17
compare_run "search_exact" \
    --search_exact "${QUERY}" --db "${DB}" --alnout out.aln --blast6out out.b6 \
    --uc out.uc --userout out.user --userfields 'query+target+id+alnlen+qrow+trow'
compare_run "allpairs_global" \
    --allpairs_global "${QUERY}" --id 0.5 --alnout out.aln --blast6out out.b6 \
    --uc out.uc --userout out.user --userfields 'query+target+id+caln'
compare_run "sintax" \
    --sintax "${QUERY}" --db "${TAXDB}" --tabbedout out.tab
compare_run "sintax --sintax_random" \
    --sintax "${QUERY}" --db "${TAXDB}" --tabbedout out.tab --sintax_random --randseed 7

# clustering: cluster.cpp, msa.cpp, otutable.cpp, linmemalign cigar
for cluster_command in --cluster_fast --cluster_size --cluster_smallmem ; do
    # --usersort applies to --cluster_smallmem only, which is the one command
    # that does not sort its input for itself
    usersort=()
    [[ "${cluster_command}" == --cluster_smallmem ]] && usersort=(--usersort)
    compare_run "${cluster_command}" \
        "${cluster_command}" "${QUERY}" --id 0.9 "${usersort[@]}" \
        --centroids out.cent --uc out.uc --msaout out.msa --consout out.cons \
        --profile out.prof --otutabout out.otu --biomout out.biom \
        --mothur_shared_out out.mothur --sizeout --clusterout_id
done
compare_run "cluster_size --clusterout_sort" \
    --cluster_size "${QUERY}" --id 0.9 --msaout out.msa --consout out.cons \
    --profile out.prof --clusterout_sort --clusterout_id --sizein --sizeout
compare_run "cluster_unoise" \
    --cluster_unoise "${QUERY}" --minsize 1 --sizein --centroids out.cent \
    --uc out.uc --otutabout out.otu
compare_run "cluster_fast --strand both --iddef 0" \
    --cluster_fast "${QUERY}" --id 0.8 --strand both --iddef 0 \
    --uc out.uc --centroids out.cent --sizeout

# chimeras: chimera.cpp (uchime and the new algorithm)
for chimera_command in --uchime_denovo --uchime2_denovo --uchime3_denovo ; do
    compare_run "${chimera_command}" \
        "${chimera_command}" "${CHIM}" --uchimeout out.uchimeout \
        --uchimealns out.alns --chimeras out.chim --nonchimeras out.nonchim \
        --borderline out.borderline --sizein --sizeout --abskew 1.5
done
compare_run "uchime_denovo --uchimeout5" \
    --uchime_denovo "${CHIM}" --uchimeout out.uchimeout --uchimeout5 \
    --fasta_score --chimeras out.chim --abskew 1.5
compare_run "uchime_ref" \
    --uchime_ref "${CHIM}" --db "${DB}" --uchimeout out.uchimeout \
    --uchimealns out.alns --nonchimeras out.nonchim
compare_run "chimeras_denovo" \
    --chimeras_denovo "${CHIM}" --tabbedout out.tab --alnout out.aln \
    --chimeras out.chim --nonchimeras out.nonchim --sizein --sizeout

# dereplication: derep.cpp, derep_prefix, derep_smallmem, fasta.cpp
compare_run "derep_fulllength" \
    --derep_fulllength "${QUERY}" --output out.fas --uc out.uc \
    --sizein --sizeout --minuniquesize 1 --strand both
compare_run "derep_fulllength --relabel" \
    --derep_fulllength "${QUERY}" --output out.fas --uc out.uc --sizeout \
    --relabel 'seq' --relabel_keep
compare_run "derep_fulllength --relabel_md5" \
    --derep_fulllength "${QUERY}" --output out.fas --relabel_md5 --sizeout
compare_run "derep_fulllength --relabel_sha1" \
    --derep_fulllength "${QUERY}" --output out.fas --relabel_sha1 --sizeout
compare_run "derep_fulllength --relabel_self" \
    --derep_fulllength "${QUERY}" --output out.fas --relabel_self --sizeout
compare_run "derep_prefix" \
    --derep_prefix "${PREFIX}" --output out.fas --uc out.uc --sizein --sizeout
compare_run "derep_smallmem" \
    --derep_smallmem "${QUERY}" --fastaout out.fas --sizein --sizeout
compare_run "fastx_uniques" \
    --fastx_uniques "${FASTQ}" --fastqout out.fastq --fastaout out.fas \
    --uc out.uc --sizein --sizeout
compare_run "rereplicate" \
    --rereplicate "${QUERY}" --output out.fas --sizein --sizeout

# FASTA/FASTQ record writers and the statistics blocks
compare_run "fastq_chars" --fastq_chars "${FASTQ}"
# --fastq_stats writes its table to the log and to stderr, not to a file
compare_run "fastq_stats" --fastq_stats "${FASTQ}"
compare_run "fastq_eestats" --fastq_eestats "${FASTQ}" --output out.eestats
compare_run "fastq_eestats2" --fastq_eestats2 "${FASTQ}" --output out.eestats
compare_run "fastq_eestats2 --length_cutoffs" \
    --fastq_eestats2 "${FASTQ}" --output out.eestats --length_cutoffs 10,80,10
compare_run "fastq_convert" \
    --fastq_convert "${FASTQ}" --fastqout out.fastq --fastq_ascii 33 --fastq_asciiout 64
compare_run "fasta2fastq" --fasta2fastq "${QUERY}" --fastqout out.fastq
compare_run "fastq_filter" \
    --fastq_filter "${FASTQ}" --fastqout out.fastq --fastaout out.fas \
    --fastq_maxee 1.0 --fastq_trunclen 40 --fastq_qmax 41 --sizeout --relabel 'f'
compare_run "fastx_filter" \
    --fastx_filter "${FASTQ}" --fastqout out.fastq --fastaout out.fas \
    --fastq_stripleft 3 --fastq_stripright 3 --fastq_qmax 41
compare_run "fastq_join" \
    --fastq_join "${R1}" --reverse "${R2}" --fastqout out.fastq \
    --join_padgap NNNN --join_padgapq 'IIII'
compare_run "fastq_mergepairs" \
    --fastq_mergepairs "${R1}" --reverse "${R2}" --fastqout out.merged \
    --fastaout out.fas --fastqout_notmerged_fwd out.notmerged1 \
    --fastqout_notmerged_rev out.notmerged2 --eetabbedout out.tab \
    --fastq_allowmergestagger
compare_run "fastx_syncpairs" \
    --fastx_syncpairs "${R1}" --reverse "${R2}" --fastqout out.r1 --fastqout_rev out.r2
compare_run "fastx_subsample" \
    --fastx_subsample "${QUERY}" --sample_pct 50 --randseed 3 --fastaout out.fas \
    --fastaout_discarded out.rejected --sizein --sizeout

# whole-file transforms and the small commands
compare_run "cut" --cut "${QUERY}" --cut_pattern 'CC^CGG_G' --fastaout out.fas \
    --fastaout_rev out.r1 --fastaout_discarded out.rejected
compare_run "fastx_mask" --fastx_mask "${QUERY}" --fastaout out.fas --hardmask
compare_run "fastx_mask --qmask none" \
    --fastx_mask "${QUERY}" --fastaout out.fas --qmask none
compare_run "maskfasta" --maskfasta "${QUERY}" --output out.fas
compare_run "fastx_revcomp" --fastx_revcomp "${QUERY}" --fastaout out.fas --label_suffix '_rc'
compare_run "sortbylength" --sortbylength "${QUERY}" --output out.fas --sizeout
compare_run "sortbysize" --sortbysize "${QUERY}" --output out.fas --sizein --sizeout
compare_run "shuffle" --shuffle "${QUERY}" --output out.fas --randseed 11
compare_run "fastx_getseq" --fastx_getseq "${DB}" --label ref1 --fastaout out.fas \
    --notmatched out.notmatched
compare_run "fastx_getseqs" --fastx_getseqs "${DB}" --labels "${LABELS}" \
    --fastaout out.fas --notmatched out.notmatched
compare_run "fastx_getsubseq" --fastx_getsubseq "${DB}" --label ref1 \
    --subseq_start 5 --subseq_end 30 --fastaout out.fas
compare_run "fasta_width 0" --derep_fulllength "${QUERY}" --output out.fas --fasta_width 0
compare_run "fasta_width 7" --derep_fulllength "${QUERY}" --output out.fas --fasta_width 7

# databases and the informational commands
compare_run "makeudb_usearch" --makeudb_usearch "${DB}" --output out.txt
compare_run "udbinfo" --udbinfo "${UDB}"
compare_run "udbstats" --udbstats "${UDB}"
compare_run "udb2fasta" --udb2fasta "${UDB}" --output out.fas
compare_run "version" --version
compare_run "help" --help
compare_run "no command" --quiet

# error paths: cli.cc's diagnostics and fatal()
compare_run "unknown option" --this_option_does_not_exist
compare_run "missing input" --derep_fulllength "${WORKDIR}/absent.fas" --output out.fas
compare_run "bad id value" --usearch_global "${QUERY}" --db "${DB}" --id 7
compare_run "empty input" --derep_fulllength /dev/null --output out.fas
compare_run "fastq of a fasta" --fastq_chars "${QUERY}"

echo "---"
echo "${comparisons} comparisons, ${differences} differences (--threads ${THREADS})"
# Runs that fail agree with each other just as well as runs that succeed, so
# they are named rather than only counted: a matrix where a whole family
# aborted is reported as identical and is worth nothing there. Only the five
# error-path entries at the end of the matrix are expected in this list.
echo "${failed_runs} of ${runs} vsearch runs exited non-zero:"
printf '  %s\n' "${FAILED_LABELS[@]}" | sort -u
[[ "${differences}" -eq 0 ]] && echo "BYTE-IDENTICAL" || echo "NOT IDENTICAL"
exit "${differences}"
