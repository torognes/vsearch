% vsearch-orient(1) version 2.30.6 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./fragments/date.md)

# NAME

vsearch \-\-orient --- use a reference database to orient fasta or fastq sequences


# SYNOPSIS

| **vsearch** **\-\-orient** _fastxfile_ **\-\-db** (_fastxfile_ | _udbfile_) (**\-\-fastaout** | **\-\-fastqout** | **\-\-notmatched** | **\-\-tabbedout**) _outputfile_ \[_options_]


# DESCRIPTION

The vsearch command `--orient` detects the orientation of input sequences by
comparing them to a reference database specified with `--db`. The two strands
of each input sequence are decomposed into words (*k*-mers) of length
`--wordlength` (default 12) and compared against words from the database. If
one strand matches at least 4 times as many words as the other, that strand is
selected as the correct orientation. Otherwise, the orientation is reported as
undetermined (?).

Correctly oriented sequences are written with `--fastaout` or `--fastqout`.
Sequences with undetermined orientation are written with `--notmatched`. A
summary table can be written with `--tabbedout`. At least one output option is
required.

To illustrate the three possible outcomes for an input sequence:

```text
Outcome  Symbol  Meaning
-------  ------  --------------------------------------------------
forward  +       forward strand matches the database
reverse  -       reverse-complementary strand matches the database
unknown  ?       neither strand clearly matches (written to --notmatched)
```


# OPTIONS

## mandatory options

#(./fragments/option_db_orient.md)

At least one of the following output options is required:

#(./fragments/option_fastaout_orient.md)

#(./fragments/option_fastqout_orient.md)

#(./fragments/option_notmatched_orient.md)

#(./fragments/option_tabbedout_orient.md)


## core options

#(./fragments/option_wordlength.md)


## secondary options

#(./fragments/option_bzip2_decompress.md)

#(./fragments/option_dbmask.md)

#(./fragments/option_fasta_width.md)

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

Use a reference database to orient sequences from a fastq file:

```sh
vsearch \
    --orient query.fastq \
    --db db.fastq \
    --fastaout query_oriented.fasta
```

Search for mis-oriented sequences by comparing a database against itself,
writing a summary table:

```sh
vsearch \
    --orient db.fasta \
    --db db.fasta \
    --tabbedout db_outliers.tsv
```

Orient sequences and collect those with undetermined orientation separately:

```sh
vsearch \
    --orient query.fasta \
    --db db.fasta \
    --fastaout oriented.fasta \
    --notmatched undetermined.fasta
```


# SEE ALSO

[`vsearch-makeudb_usearch(1)`](./vsearch-makeudb_usearch.1.md),
[`vsearch-udb(5)`](../formats/vsearch-udb.5.md)


#(./fragments/footer.md)
