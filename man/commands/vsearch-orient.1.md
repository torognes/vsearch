% vsearch-orient(1) version 2.30.0 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./fragments/date.md)

# NAME

vsearch \-\-orient --- use a reference database to orient fasta or
fastq sequences


# SYNOPSIS

| **vsearch** **\-\-orient** _fastxfile_ \-\-db (_fastxfile_ | _udbfile_) (\-\-fastaout | \-\-fastqout | \-\-notmatched | \-\-tabbedout) _outputfile_ \[_options_]


# DESCRIPTION

The vsearch command `--orient` detects if input sequences are in the
same orientation as sequences in the reference database specified with
the `--db` option. Detection is based on the number of common *words*
(*k*-mers) of length 12.

Each reference sequence is decomposed into a set of words. The forward
(+) and reverse-complement (-) strands of each input sequence are also
decomposed into sets of words. If a strand shares 4 times (or more)
words with a reference sequence than the other strand, then it is
chosen. Otherwise, the input sequence is marked as undetermined (?).


# OPTIONS

## mandatory options

`--db` *filename*
: Read the reference database from the given fasta, fastq, or UDB
  *filename*. Only UDB files created with a `--wordlength` of 12 are
  accepted (see
  [`vsearch-makeudb_usearch(1)`](./vsearch-makeudb_usearch.1.md) for
  more details).

At least one of `--fastaout`, `--fastqout`, `--notmatched`,
`--tabbedout` must also be specified.


## core options

`--fastaout` *filename*
: Write the correctly oriented sequences to *filename*, in fasta
  format.

`--fastqout` *filename*
: Write the correctly oriented sequences to *filename*, in fastq
  format (requires a fastq input).

`--notmatched` *filename*
: Write the sequences with undetermined direction to *filename*, in
  their original format.

`--tabbedout` *filename*
: Write the resuls to a four-column tab-delimited file with the
  specified *filename*. Columns are:

    1. query label
    2. direction (+, -, or ?)
    3. number of matching words on the forward strand
    4. number of matching words on the reverse complementary strand


## secondary options

#(./fragments/option_bzip2_decompress.md)

#(./fragments/option_fasta_width.md)

#(./fragments/option_dbmask.md)

#(./fragments/option_gzip_decompress.md)

#(./fragments/option_label_suffix.md)

#(./fragments/option_lengthout.md)

#(./fragments/option_log.md)

#(./fragments/option_no_progress.md)

#(./fragments/option_notrunclabels.md)

#(./fragments/option_qmask.md)

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

Use the sequences in *db.fastq* to orient the sequences in
*query.fastq*. Write the results to *query_oriented.fasta*, in fasta
format:

```sh
vsearch \
    --orient query.fastq \
    --db db.fastq \
    --fastaout query_oriented.fasta
```

Search for mis-oriented in a set of reference sequences *db.fasta* by
comparing it to itself. Write the results to the table
*db_outliers.table*:

```sh
vsearch \
    --orient db.fasta \
    --db db.fasta \
    --tabbedout db_outliers.table
```


# SEE ALSO

[`vsearch-makeudb_usearch(1)`](./vsearch-makeudb_usearch.1.md),
[`vsearch-udb(5)`](./vsearch-udb.5.md)


#(./fragments/footer.md)
