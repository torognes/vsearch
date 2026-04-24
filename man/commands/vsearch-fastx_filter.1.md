% vsearch-fastx_filter(1) version 2.31.0 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./fragments/date.md)

# NAME

vsearch \-\-fastx_filter --- trim and filter fasta or fastq sequences


# SYNOPSIS

| **vsearch** **\-\-fastx_filter** _inputfile_ \[**\-\-reverse** _inputfile_] (**\-\-fastaout** | **\-\-fastaout_discarded** | **\-\-fastaout_discarded_rev** | **\-\-fastaout_rev** | **\-\-fastqout** | **\-\-fastqout_discarded** | **\-\-fastqout_discarded_rev** | **\-\-fastqout_rev**) _filename_ \[_options_]


# DESCRIPTION

The vsearch command `--fastx_filter` trims and filters the sequences
in a fasta or fastq file. The input format is detected automatically
(see [`vsearch-fasta(5)`](../formats/vsearch-fasta.5.md) and
[`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md)). Sequences that
pass all filters are written to the output files specified with
`--fastaout` and/or `--fastqout`. Discarded sequences can be written
to `--fastaout_discarded` and/or `--fastqout_discarded`. Output
cannot be written to fastq files when the input is in fasta format,
since quality scores are not available.

Sequences are first **trimmed**, then **filtered** based on the
remaining bases:

- **Trimming** shortens sequences using `--fastq_stripleft`,
  `--fastq_stripright`, `--fastq_truncee`, `--fastq_truncee_rate`,
  `--fastq_trunclen`, `--fastq_trunclen_keep`, and
  `--fastq_truncqual`.

- **Filtering** discards sequences that do not satisfy criteria set by
  `--fastq_maxee`, `--fastq_maxee_rate`, `--fastq_maxlen`,
  `--fastq_maxns`, `--fastq_minlen` (default 1), `--fastq_minqual`,
  `--fastq_trunclen`, `--maxsize`, and `--minsize`.

If no trimming or filtering options are given, all sequences are
written to the output files, possibly after conversion from fastq to
fasta format.

For paired-end reads, the file with reverse reads is specified with
`--reverse`, and the corresponding outputs are written to
`--fastaout_rev`, `--fastqout_rev`, `--fastaout_discarded_rev`, and
`--fastqout_discarded_rev`. Both reads of a pair must pass all filters
for either to be retained; if one fails, both are discarded.

When the input is in fasta format, the following options are not
accepted because quality scores are not available: `--eeout`,
`--fastq_eeout`, `--fastq_maxee`, `--fastq_maxee_rate`,
`--fastq_minqual`, `--fastq_truncee`, `--fastq_truncee_rate`,
`--fastqout`, `--fastqout_discarded`, `--fastqout_discarded_rev`, and
`--fastqout_rev`. The following options are silently ignored when
reading fasta: `--fastq_ascii`, `--fastq_qmax`, `--fastq_qmin`, and
`--fastq_truncqual`.

Sequence headers are not truncated at whitespace by this command,
regardless of input format. The `--notrunclabels` option is therefore
effectively a no-op with `--fastx_filter`.

After processing, vsearch reports the number of sequences kept and
discarded, and how many of the kept sequences were trimmed. Use
`--eeout` or `--fastq_eeout` to annotate output headers with the
expected error count. See
[`vsearch-expected_error(7)`](../misc/vsearch-expected_error.7.md)
for details on expected errors.

The `--fastq_filter` command is similar but restricted to fastq input
(see [`vsearch-fastq_filter(1)`](./vsearch-fastq_filter.1.md)).

To illustrate the effect of `--fastq_minlen 5` on a fasta input:

```text
Input (3 reads):   --fastq_minlen 5:

>r1                kept (--fastaout):      discarded (--fastaout_discarded):
ACGTACGT           >r1  ACGTACGT          >r3 (4 bases < 5)
>r2         -->    >r2  ACGTAC
ACGTAC
>r3
ACGT
```


# OPTIONS

## mandatory options

At least one output option must be specified: `--fastaout`,
`--fastaout_discarded`, `--fastaout_discarded_rev`, `--fastaout_rev`,
`--fastqout`, `--fastqout_discarded`, `--fastqout_discarded_rev`, or
`--fastqout_rev`.


## core options

### trimming

#(./fragments/option_fastq_stripleft.md)

#(./fragments/option_fastq_stripright.md)

#(./fragments/option_fastq_truncee.md)

#(./fragments/option_fastq_truncee_rate.md)

#(./fragments/option_fastq_trunclen.md)

#(./fragments/option_fastq_trunclen_keep.md)

#(./fragments/option_fastq_truncqual.md)


### filtering

#(./fragments/option_fastq_maxee.md)

#(./fragments/option_fastq_maxee_rate.md)

#(./fragments/option_fastq_maxlen.md)

#(./fragments/option_fastq_maxns.md)

#(./fragments/option_fastq_minlen.md)

#(./fragments/option_fastq_minqual.md)

#(./fragments/option_maxsize.md)

#(./fragments/option_minsize.md)


### output

#(./fragments/option_eeout.md)

#(./fragments/option_fastq_eeout.md)

#(./fragments/option_fastaout.md)

#(./fragments/option_fastaout_discarded.md)

#(./fragments/option_fastaout_discarded_rev.md)

#(./fragments/option_fastaout_rev.md)

#(./fragments/option_fastqout.md)

#(./fragments/option_fastqout_discarded.md)

#(./fragments/option_fastqout_discarded_rev.md)

#(./fragments/option_fastqout_rev.md)

#(./fragments/option_reverse.md)

#(./fragments/option_sizein.md)

#(./fragments/option_sizeout.md)


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

#(./fragments/option_xee.md)

#(./fragments/option_xlength.md)

#(./fragments/option_xsize.md)


## ignored options

#(./fragments/option_threads_not_multithreaded.md)


# EXAMPLES

Filter a fastq file, discarding reads with an expected error above
1.0, and write both kept and discarded sequences:

```sh
vsearch \
    --fastx_filter input.fastq \
    --fastq_maxee 1.0 \
    --fastqout kept.fastq \
    --fastqout_discarded discarded.fastq
```

Truncate reads to 200 bases, then discard those with an expected error
above 0.5, and annotate output headers with the expected error count:

```sh
vsearch \
    --fastx_filter input.fastq \
    --fastq_trunclen 200 \
    --fastq_maxee 0.5 \
    --eeout \
    --fastqout filtered.fastq
```

Strip primers from both ends of reads in a fasta file, discarding
reads shorter than 50 bases after stripping:

```sh
vsearch \
    --fastx_filter input.fasta \
    --fastq_stripleft 20 \
    --fastq_stripright 20 \
    --fastq_minlen 50 \
    --fastaout stripped.fasta
```

Filter paired-end reads, keeping only pairs where both reads pass:

```sh
vsearch \
    --fastx_filter forward.fastq \
    --reverse reverse.fastq \
    --fastq_maxee 1.0 \
    --fastqout kept_fwd.fastq \
    --fastqout_rev kept_rev.fastq \
    --fastqout_discarded discarded_fwd.fastq \
    --fastqout_discarded_rev discarded_rev.fastq
```


# SEE ALSO

[`vsearch-fastq_filter(1)`](./vsearch-fastq_filter.1.md),
[`vsearch-fastq_chars(1)`](./vsearch-fastq_chars.1.md),
[`vsearch-fastq_stats(1)`](./vsearch-fastq_stats.1.md),
[`vsearch-fasta(5)`](../formats/vsearch-fasta.5.md),
[`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md),
[`vsearch-expected_error(7)`](../misc/vsearch-expected_error.7.md)


#(./fragments/footer.md)
