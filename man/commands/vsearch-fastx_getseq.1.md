% vsearch-fastx_getseq(1) version 2.30.4 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./fragments/date.md)

# NAME

vsearch \-\-fastx_getseq --- extract a sequence from a fasta or fastq
file by label


# SYNOPSIS

| **vsearch** **\-\-fastx_getseq** _fastxfile_ (**\-\-fastaout** | **\-\-fastqout** | **\-\-notmatched** | **\-\-notmatchedfq**) _filename_ **\-\-label** _string_ \[_options_]


# DESCRIPTION

The vsearch command `--fastx_getseq` extracts a sequence from a fasta
or fastq file whose header matches the label given with `--label`. The
input format is detected automatically (see
[`vsearch-fasta(5)`](../formats/vsearch-fasta.5.md) and
[`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md)).

By default, the label must match the entire header (after truncation
at the first space or tab — use `--notrunclabels` to disable
truncation). Use `--label_substr_match` to allow the label to match
anywhere in the header instead. Matching is not case-sensitive.

Matching sequences are written to `--fastaout` and/or `--fastqout`.
Non-matching sequences can be written to `--notmatched` and/or
`--notmatchedfq`.

To extract multiple sequences, or to match by word rather than full
label, see
[`vsearch-fastx_getseqs(1)`](./vsearch-fastx_getseqs.1.md). To
extract a subsequence region, see
[`vsearch-fastx_getsubseq(1)`](./vsearch-fastx_getsubseq.1.md).

```text
Input (3 sequences):  --label "seq2":

>seq1                 extracted (--fastaout):   not matched (--notmatched):
AAAA                  >seq2  CCCC              >seq1  AAAA
>seq2         -->                               >seq3  TTTT
CCCC
>seq3
TTTT
```


# OPTIONS

## mandatory options

`--label` is required. At least one of `--fastaout`, `--fastqout`,
`--notmatched`, or `--notmatchedfq` must be specified.


## core options

#(./fragments/option_label.md)

#(./fragments/option_label_substr_match.md)

`--fastaout` *filename*
: Write the matching sequences to *filename*, in fasta format.

`--fastqout` *filename*
: Write the matching sequences to *filename*, in fastq format (see
  [`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md)). Requires fastq
  input.

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

Extract the sequence labelled `seq42` from a fasta file:

```sh
vsearch \
    --fastx_getseq input.fasta \
    --label "seq42" \
    --fastaout extracted.fasta
```

Extract a sequence by substring match (case-insensitive), writing
non-matching sequences to a separate file:

```sh
vsearch \
    --fastx_getseq input.fastq \
    --label "primer" \
    --label_substr_match \
    --fastqout matched.fastq \
    --notmatchedfq unmatched.fastq
```


# SEE ALSO

[`vsearch-fastx_getseqs(1)`](./vsearch-fastx_getseqs.1.md),
[`vsearch-fastx_getsubseq(1)`](./vsearch-fastx_getsubseq.1.md),
[`vsearch-fasta(5)`](../formats/vsearch-fasta.5.md),
[`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md)


#(./fragments/footer.md)
