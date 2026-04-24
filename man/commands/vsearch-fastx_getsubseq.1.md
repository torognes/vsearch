% vsearch-fastx_getsubseq(1) version 2.31.0 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./fragments/date.md)

# NAME

vsearch \-\-fastx_getsubseq --- extract a subsequence from a fasta or
fastq file by label


# SYNOPSIS

| **vsearch** **\-\-fastx_getsubseq** _fastxfile_ (**\-\-fastaout** | **\-\-fastqout** | **\-\-notmatched** | **\-\-notmatchedfq**) _filename_ **\-\-label** _string_ \[**\-\-subseq_start** _position_] \[**\-\-subseq_end** _position_] \[_options_]


# DESCRIPTION

The vsearch command `--fastx_getsubseq` extracts a subsequence region
from a sequence in a fasta or fastq file whose header matches the
label given with `--label`. The input format is detected automatically
(see [`vsearch-fasta(5)`](../formats/vsearch-fasta.5.md) and
[`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md)).

This command is similar to `--fastx_getseq` (see
[`vsearch-fastx_getseq(1)`](./vsearch-fastx_getseq.1.md)), but
instead of writing the full matching sequence, it writes only the
region specified by `--subseq_start` and `--subseq_end`. Positions
are 1-based: the first base is at position 1. By default, extraction
starts at position 1 and ends at the last position of the sequence.

Label matching follows the same rules as `--fastx_getseq`: the label
must match the entire header by default (after truncation at the first
space or tab — use `--notrunclabels` to disable truncation). Use
`--label_substr_match` to allow matching anywhere in the header.
Matching is not case-sensitive.

Matching sequences (trimmed to the requested region) are written to
`--fastaout` and/or `--fastqout`. Non-matching sequences are written
in full to `--notmatched` and/or `--notmatchedfq`.

```text
Input (2 sequences):  --label "seq1" --subseq_start 3 --subseq_end 6:

>seq1                 extracted (--fastaout):   not matched (--notmatched):
ACGTACGT      -->     >seq1  GTAC              >seq2  TTTTTTTT
>seq2
TTTTTTTT
```


# OPTIONS

## mandatory options

`--label` is required. At least one of `--fastaout`, `--fastqout`,
`--notmatched`, or `--notmatchedfq` must be specified.


## core options

#(./fragments/option_label.md)

#(./fragments/option_label_substr_match.md)

#(./fragments/option_subseq_start.md)

#(./fragments/option_subseq_end.md)

`--fastaout` *filename*
: Write the extracted subsequences to *filename*, in fasta format.

`--fastqout` *filename*
: Write the extracted subsequences to *filename*, in fastq format
  (see [`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md)). Requires
  fastq input.

#(./fragments/option_notmatched.md)

#(./fragments/option_notmatchedfq.md)


## secondary options

#(./fragments/option_bzip2_decompress.md)

#(./fragments/option_fasta_width.md)

#(./fragments/option_fastq_ascii.md)

#(./fragments/option_fastq_qmax.md)

#(./fragments/option_fastq_qmin.md)

#(./fragments/option_gzip_decompress.md)

#(./fragments/option_label_suffix.md)

#(./fragments/option_lengthout.md)

#(./fragments/option_log.md)

#(./fragments/option_no_progress.md)

#(./fragments/option_notrunclabels.md)

#(./fragments/option_quiet.md)

#(./fragments/option_relabel.md)

#(./fragments/option_relabel_keep.md)

#(./fragments/option_relabel_md5.md)

#(./fragments/option_relabel_self.md)

#(./fragments/option_relabel_sha1.md)

#(./fragments/option_sample.md)

#(./fragments/option_sizein.md)

#(./fragments/option_sizeout.md)

#(./fragments/option_xee.md)

#(./fragments/option_xlength.md)

#(./fragments/option_xsize.md)


## ignored options

#(./fragments/option_threads_not_multithreaded.md)


# EXAMPLES

Extract bases 21 to 270 from a target sequence, discarding the rest:

```sh
vsearch \
    --fastx_getsubseq input.fasta \
    --label "target_seq" \
    --subseq_start 21 \
    --subseq_end 270 \
    --fastaout region.fasta
```

Extract from position 100 to the end of the sequence in a fastq file:

```sh
vsearch \
    --fastx_getsubseq input.fastq \
    --label "read1" \
    --subseq_start 100 \
    --fastqout tail.fastq \
    --notmatchedfq unmatched.fastq
```


# SEE ALSO

[`vsearch-fastx_getseq(1)`](./vsearch-fastx_getseq.1.md),
[`vsearch-fastx_getseqs(1)`](./vsearch-fastx_getseqs.1.md),
[`vsearch-fasta(5)`](../formats/vsearch-fasta.5.md),
[`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md)


#(./fragments/footer.md)
