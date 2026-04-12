% vsearch-fastq_filter(1) version 2.30.6 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./fragments/date.md)

# NAME

vsearch \-\-fastq_filter --- trim and filter fastq sequences


# SYNOPSIS

| **vsearch** **\-\-fastq_filter** _fastqfile_ \[**\-\-reverse** _fastqfile_] (**\-\-fastaout** | **\-\-fastaout_discarded** | **\-\-fastaout_discarded_rev** | **\-\-fastaout_rev** | **\-\-fastqout** | **\-\-fastqout_discarded** | **\-\-fastqout_discarded_rev** | **\-\-fastqout_rev**) _filename_ \[_options_]


# DESCRIPTION

The vsearch command `--fastq_filter` trims and filters the sequences
in a fastq file (see
[`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md)). It is similar
to `--fastx_filter` (see
[`vsearch-fastx_filter(1)`](./vsearch-fastx_filter.1.md)), but
restricted to fastq input only.

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
written to the output files, possibly after conversion to fasta
format.

For paired-end reads, the file with reverse reads is specified with
`--reverse`, and the corresponding outputs are written to
`--fastaout_rev`, `--fastqout_rev`, `--fastaout_discarded_rev`, and
`--fastqout_discarded_rev`. Both reads of a pair must pass all filters
for either to be retained; if one fails, both are discarded.

After processing, vsearch reports the number of sequences kept and
discarded, and how many of the kept sequences were trimmed. Use
`--eeout` or `--fastq_eeout` to annotate output headers with the
expected error count. See
[`vsearch-expected_error(7)`](../misc/vsearch-expected_error.7.md)
for details on expected errors.

To illustrate the effect of `--fastq_trunclen 6 --fastq_maxee 1.0`:

```text
Input (3 reads):   --fastq_trunclen 6 --fastq_maxee 1.0:

@r1                kept (--fastqout):
ACGTACGTAC         @r1 (truncated to ACGTAC; expected error < 1.0)
@r2         -->
ACGT               discarded (--fastqout_discarded):
@r3                @r2 (4 bases < 6, discarded by --fastq_trunclen)
NNNNNNNNNN         @r3 (expected error > 1.0 after truncation)
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
1.0, and write kept sequences in fasta format:

```sh
vsearch \
    --fastq_filter input.fastq \
    --fastq_maxee 1.0 \
    --fastaout filtered.fasta
```

Truncate reads to 250 bases, then discard those with an expected error
above 0.5, annotating output headers with the expected error:

```sh
vsearch \
    --fastq_filter input.fastq \
    --fastq_trunclen 250 \
    --fastq_maxee 0.5 \
    --eeout \
    --fastqout filtered.fastq
```

Strip a fixed-length primer from the left, then filter by minimum
length and maximum expected error rate:

```sh
vsearch \
    --fastq_filter input.fastq \
    --fastq_stripleft 20 \
    --fastq_minlen 100 \
    --fastq_maxee_rate 0.01 \
    --fastqout trimmed.fastq \
    --fastqout_discarded discarded.fastq
```

Filter paired-end reads, keeping only pairs where both reads pass:

```sh
vsearch \
    --fastq_filter forward.fastq \
    --reverse reverse.fastq \
    --fastq_maxee 1.0 \
    --fastqout kept_fwd.fastq \
    --fastqout_rev kept_rev.fastq \
    --fastqout_discarded discarded_fwd.fastq \
    --fastqout_discarded_rev discarded_rev.fastq
```


# SEE ALSO

[`vsearch-fastx_filter(1)`](./vsearch-fastx_filter.1.md),
[`vsearch-fastq_chars(1)`](./vsearch-fastq_chars.1.md),
[`vsearch-fastq_stats(1)`](./vsearch-fastq_stats.1.md),
[`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md),
[`vsearch-expected_error(7)`](../misc/vsearch-expected_error.7.md)


#(./fragments/footer.md)
