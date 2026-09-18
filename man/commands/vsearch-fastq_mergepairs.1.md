% vsearch-fastq_mergepairs(1) version 2.32.0 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./fragments/date.md)

# NAME

vsearch \-\-fastq_mergepairs --- merge paired-end reads into one sequence


# SYNOPSIS

| **vsearch** **\-\-fastq_mergepairs** _fwdfile_ **\-\-reverse** _revfile_ (**\-\-fastaout** | **\-\-fastqout** | _other output options_) _filename_ \[_options_]


# DESCRIPTION

The vsearch command `--fastq_mergepairs` merges paired-end sequence reads
into a single sequence by aligning the forward and reverse reads and combining
their overlapping regions. The forward reads are specified as the argument to
this option; the reverse reads are specified with `--reverse`. Reads are
matched by position: the first forward read is paired with the first reverse
read, the second with the second, and so on. Labels are not used for
matching; the run stops with a fatal error if the two files contain
different numbers of reads.

The reverse read is reverse-complemented before alignment. Merging requires
an overlap between the two reads of at least `--fastq_minovlen` bases
(default 10, minimum 5). Read pairs with too many mismatches in the overlap
— more than `--fastq_maxdiffs` (default 10) or more than
`--fastq_maxdiffpct` percent (default 100.0%) — are discarded. Additional
heuristics prevent merging of read pairs that cannot be aligned reliably.

In the merged region, quality scores from the two reads are combined using
the Phred score formula. Outside the overlap, the quality scores from the
contributing read are used directly. The quality score of `N` bases is
replaced by the minimum score (Q0) in every output, including the
not-merged output files. Output quality scores can be clamped
with `--fastq_qmaxout` and `--fastq_qminout` (these apply only to the
merged region). Unlike the commands that pass an input quality through,
`--fastq_mergepairs` keeps the pre-2.32.0 `--fastq_qmaxout` default of 41,
because the score it clamps is computed rather than read: two agreeing
Q40 bases have a posterior quality of Q85, and reporting it would change
the merged output of every run. Pass `--fastq_qmaxout 93` for the
unclamped posterior. The merged scores are written with the *input*
offset, `--fastq_ascii`: this is the only command whose output encoding
tracks its input encoding, and it does not accept `--fastq_asciiout`.

Staggered pairs — where the 3' end of one read extends past the 5' end of
the other — are discarded by default. Use `--fastq_allowmergestagger` to
allow them; the overhanging portions of both reads are excluded from the
merged sequence, which then covers the overlap only.

Reads can be pre-filtered with `--fastq_truncqual`, `--fastq_minlen`,
`--fastq_maxlen` (the length bounds apply after truncation), and
`--fastq_maxns`; the expected-error filter `--fastq_maxee` applies to
the *merged* sequence. Bounds on the merged sequence length are set
with `--fastq_minmergelen` and `--fastq_maxmergelen`.

To illustrate a merge with a 6-base overlap:

```text
Forward (5'→3'):   AAAATTTTTT
Reverse (5'→3'):   GGGGAAAAAA

Aligned:
  Forward:      AAAATTTTTT
  Rev-comp:         TTTTTTCCCC   (reverse complement of the reverse read)

Overlap:            TTTTTT      (6 bases; scores combined)
Merged:         AAAATTTTTTCCCC
```

At the end of the run, vsearch prints a report — to standard error
(unless `--quiet` is given), or to the log file when `--log` is used —
giving the number of merged
pairs and, for the pairs that could not be merged, a breakdown by
reason. Most reasons correspond directly to a user-adjustable
threshold: `reads too short (after truncation)` (`--fastq_minlen`),
`reads too long (after truncation)` (`--fastq_maxlen`), `too many
N's` (`--fastq_maxns`), `too many differences` (`--fastq_maxdiffs`),
`too high percentage of differences` (`--fastq_maxdiffpct`), `overlap
too short` (`--fastq_minovlen`), `expected error too high`
(`--fastq_maxee`), `merged fragment too short` (`--fastq_minmergelen`),
`merged fragment too long` (`--fastq_maxmergelen`), and `staggered
read pairs` (allowed with `--fastq_allowmergestagger`). The remaining
reasons reflect the alignment heuristics and are not directly
controllable:

`too few kmers found on same diagonal`
: Too few k-mers shared by the two reads fall on a common alignment
  diagonal, so no candidate overlap could be located and no alignment
  was attempted.

`multiple potential alignments`
: More than one overlap of comparable quality was found — for instance
  when the overlap region contains a tandem repeat — so the correct
  alignment is ambiguous and the pair is left unmerged.

`alignment score too low, or score drop too high`
: An overlap was found and aligned, but its score remained below the
  internal threshold, or the score dropped too sharply within the
  overlap. This usually points to clustered mismatches or an indel in
  the overlap region.


# ALGORITHM

Merging happens in three stages: candidate overlaps are located with
shared words, each candidate is scored by an ungapped alignment whose
scores are derived from the base qualities, and the best candidate is
accepted only if it is both good enough and unambiguous. The thresholds
below are internal and have no options of their own.

**Locating candidates.** The forward read is indexed by its 5-mers (a
shorter word than the search commands use), and every possible overlap
length is rated by the number of 5-mers the two reads share at that
offset. Only offsets sharing at least four 5-mers are examined further.
A pair where no offset reaches that count is discarded as `too few kmers
found on same diagonal`, without any alignment being attempted. When
`--fastq_minovlen` is below 9 the requirement is relaxed to
`--fastq_minovlen` minus 4.

**Scoring a candidate.** Each surviving offset is aligned without gaps,
the forward read against the reverse complement of the reverse read, and
scored in bits as a log-odds ratio. Writing *p* for the probability that
two truly identical bases are observed as a match, given the error
probabilities *p_f* and *p_r* of the forward and reverse base:

```text
p = 1 - p_f - p_r + 4 x p_f x p_r / 3
```

a matching column scores log2(*p* / 0.25) and a mismatching column
log2((1 - *p*) / 0.75), the latter capped at -4 bits so that one
disagreement between two poor bases cannot dominate. Two good bases that
agree are therefore worth almost 2 bits each, while a disagreement
between two good bases costs many. The formulas are those of Edgar &
Flyvbjerg (2015), which also give the posterior qualities written for
the merged region.

**Accepting a candidate.** The running score and its maximum so far are
tracked as the overlap is walked, and a candidate whose score falls more
than 16 bits below its own maximum is dropped: that is how a run of
clustered mismatches is rejected even when the two ends align well. A
candidate is *acceptable* once its final score reaches 16 bits, or 1.6
times `--fastq_minovlen` when that is below 9 --- about nine matching
bases of good quality, since eight Q40 matches come to 15.998 bits. If
more than one candidate is acceptable the overlap is ambiguous and the
pair is discarded as `multiple potential alignments`, which is what a
tandem repeat in the overlap produces. Otherwise the highest-scoring
candidate is kept, and a pair whose best candidate never reached the
threshold is discarded as `alignment score too low, or score drop too
high`.

Because the alignment is ungapped, an indel in the overlap is not
modelled as such: it shows up as a run of mismatches and is rejected by
the score-drop rule.


# OPTIONS

## mandatory options

#(./fragments/option_reverse.md)

At least one output option is required: one of the merged-read
outputs below, one of the not-merged outputs
(`--fastaout_notmerged_fwd`, `--fastaout_notmerged_rev`,
`--fastqout_notmerged_fwd`, `--fastqout_notmerged_rev`), or
`--eetabbedout` (see the following sections).

#(./fragments/option_fastaout_mergepairs.md)

#(./fragments/option_fastqout_mergepairs.md)


## core options

#(./fragments/option_fastq_ascii.md)

#(./fragments/option_fastq_minovlen.md)

#(./fragments/option_fastq_maxdiffs.md)

#(./fragments/option_fastq_maxdiffpct.md)

#(./fragments/option_fastq_allowmergestagger.md)

#(./fragments/option_fastq_nostagger.md)

#(./fragments/option_fastq_minmergelen.md)

#(./fragments/option_fastq_maxmergelen.md)

#(./fragments/option_fastq_qmax.md)

#(./fragments/option_fastq_qmin.md)

#(./fragments/option_fastq_qmaxout.md)

#(./fragments/option_fastq_qminout.md)


## secondary options

#(./fragments/option_bzip2_decompress.md)

#(./fragments/option_eeout.md)

#(./fragments/option_eetabbedout.md)

#(./fragments/option_fasta_width.md)

#(./fragments/option_fastaout_notmerged_fwd.md)

#(./fragments/option_fastaout_notmerged_rev.md)

#(./fragments/option_fastq_eeout.md)

#(./fragments/option_fastq_maxee.md)

#(./fragments/option_fastq_maxlen.md)

#(./fragments/option_fastq_maxns.md)

#(./fragments/option_fastq_minlen.md)

#(./fragments/option_fastq_truncqual.md)

#(./fragments/option_fastqout_notmerged_fwd.md)

#(./fragments/option_fastqout_notmerged_rev.md)

#(./fragments/option_gzip_decompress.md)

#(./fragments/option_label_suffix.md)

#(./fragments/option_lengthout.md)

#(./fragments/option_log.md)

#(./fragments/option_no_progress.md)

#(./fragments/option_quiet.md)

#(./fragments/option_relabel.md)

#(./fragments/option_relabel_keep.md)

#(./fragments/option_relabel_md5.md)

#(./fragments/option_relabel_self.md)

#(./fragments/option_relabel_sha1.md)

#(./fragments/option_sample.md)

#(./fragments/option_sizein.md)

#(./fragments/option_sizeout.md)

#(./fragments/option_threads.md)

#(./fragments/option_xee.md)

#(./fragments/option_xlength.md)

#(./fragments/option_xsize.md)


# EXAMPLES

Merge paired-end reads and write merged sequences to a fastq file:

```sh
vsearch \
    --fastq_mergepairs fwd.fastq \
    --reverse rev.fastq \
    --fastqout merged.fastq
```

Merge with a stricter overlap and mismatch threshold, and save
unmerged reads for inspection:

```sh
vsearch \
    --fastq_mergepairs fwd.fastq \
    --reverse rev.fastq \
    --fastq_minovlen 20 \
    --fastq_maxdiffs 5 \
    --fastqout merged.fastq \
    --fastqout_notmerged_fwd unmerged_fwd.fastq \
    --fastqout_notmerged_rev unmerged_rev.fastq
```

Allow staggered pairs and filter on expected error after merging:

```sh
vsearch \
    --fastq_mergepairs fwd.fastq \
    --reverse rev.fastq \
    --fastq_allowmergestagger \
    --fastq_maxee 1.0 \
    --fastqout merged.fastq
```


# SEE ALSO

[`vsearch-fastx_filter(1)`](./vsearch-fastx_filter.1.md),
[`vsearch-fastq_eestats(1)`](./vsearch-fastq_eestats.1.md),
[`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md),
[`vsearch-expected_error(7)`](../misc/vsearch-expected_error.7.md),
[`vsearch-usearch(7)`](../misc/vsearch-usearch.7.md)


#(./fragments/footer.md)
