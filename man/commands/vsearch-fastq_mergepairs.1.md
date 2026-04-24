% vsearch-fastq_mergepairs(1) version 2.31.0 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./fragments/date.md)

# NAME

vsearch \-\-fastq_mergepairs --- merge paired-end reads into one sequence


# SYNOPSIS

| **vsearch** **\-\-fastq_mergepairs** _fwdfile_ **\-\-reverse** _revfile_ (**\-\-fastaout** | **\-\-fastqout**) _filename_ \[_options_]


# DESCRIPTION

The vsearch command `--fastq_mergepairs` merges paired-end sequence reads
into a single sequence by aligning the forward and reverse reads and combining
their overlapping regions. The forward reads are specified as the argument to
this option; the reverse reads are specified with `--reverse`. Reads are
matched by position: the first forward read is paired with the first reverse
read, the second with the second, and so on. Labels are not used for
matching, but a warning is emitted if the two files contain different numbers
of reads.

The reverse read is reverse-complemented before alignment. Merging requires
an overlap between the two reads of at least `--fastq_minovlen` bases
(default 10, minimum 5). Read pairs with too many mismatches in the overlap
— more than `--fastq_maxdiffs` (default 10) or more than
`--fastq_maxdiffpct` percent (default 100.0%) — are discarded. Additional
heuristics prevent merging of read pairs that cannot be aligned reliably.

In the merged region, quality scores from the two reads are combined using
the Phred score formula. Outside the overlap, the quality scores from the
contributing read are used directly. Output quality scores can be clamped
with `--fastq_qmaxout` and `--fastq_qminout` (these apply only to the
merged region).

Staggered pairs — where the 3' end of the reverse read extends past the 5'
end of the forward read — are discarded by default. Use
`--fastq_allowmergestagger` to allow them; the overhanging portion is
excluded from the merged sequence.

Reads can be pre-filtered with `--fastq_truncqual`, `--fastq_maxee`,
`--fastq_minlen`, `--fastq_maxlen`, and `--fastq_maxns`. Bounds on the
merged sequence length are set with `--fastq_minmergelen` and
`--fastq_maxmergelen`.

To illustrate a merge with a 6-base overlap:

```text
Forward (5'→3'):   AAAAAAAAATTTTTTTGGG
Reverse (5'→3'):   CCCCCCCCCCGGGGGG         (reverse complement of rev read)

Aligned:
  Forward:    AAAAAAAAATTTTTTTGGG
  Rev-comp:         TTTTTTTGGGCCCCCCCCCC

Overlap:            TTTTTTTGGG  (6 bases; scores combined)
Merged:       AAAAAAAAATTTTTTTGGGCCCCCCCCCC
```


# OPTIONS

## mandatory options

#(./fragments/option_reverse.md)

At least one of the following output options is required:

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
[`vsearch-expected_error(7)`](../misc/vsearch-expected_error.7.md)


#(./fragments/footer.md)
