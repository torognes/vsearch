% vsearch-udbstats(1) version 2.31.0 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./fragments/date.md)

# NAME

vsearch \-\-udbstats --- report statistics about indexed words in a UDB database file


# SYNOPSIS

| **vsearch** **\-\-udbstats** _udbfile_ \[_options_]


# DESCRIPTION

The vsearch command `--udbstats` reads a UDB database file and writes
statistics about its indexed words (*k*-mers) to the log file: the
report is only produced when `--log` is given (without it, only the
generic loading summary appears on the standard error). The report
gives the distribution of word frequencies across the indexed
sequences in the database.

The report has four parts: a header (word length, number of index
slots, DBAccel), overall totals (database size in nucleotides, number
of indexed words, and the median and mean number of sequences per
word), a table of the eleven most frequent words giving each word's
index and sequence form, its number of matching sequences (Size) and
the first eight of them (Row), and a histogram binning the index slots
by how many sequences they match, in buckets that double in width.

Several reported fields are USEARCH header fields vsearch does not use,
and are constants rather than measurements: Word ones repeats the word
length, Spaced, Hashed, Coded and Stepped are always `No`, the Cap
column is always `0`, and Lower is always `0` (masking is not recorded
in a UDB file, so Upper and Total both equal the nucleotide count).

Note that the number of index slots is 4 raised to the power of the
word length, so the cost of this command is set by the word length
rather than by the size of the database.

See [`vsearch-udb(5)`](../formats/vsearch-udb.5.md) for a description
of the UDB file format.


# OPTIONS

## mandatory options

`--udbstats` *udbfile*
: Read the UDB database *udbfile* and report word statistics. As UDB
  files cannot be read from pipes, *udbfile* must be a seekable file path
  (see [`vsearch-udb(5)`](../formats/vsearch-udb.5.md)); a pipe is
  rejected with an error.


## secondary options

#(./fragments/option_log.md)

#(./fragments/option_no_progress.md)

#(./fragments/option_quiet.md)


## ignored options

#(./fragments/option_threads_not_multithreaded.md)


# EXAMPLES

Report word statistics for a UDB database (the report needs `--log`;
without it the command prints only the loading summary):

```sh
vsearch \
    --udbstats db.udb \
    --log stats.log
```


# SEE ALSO

[`vsearch-makeudb_usearch(1)`](./vsearch-makeudb_usearch.1.md),
[`vsearch-udb2fasta(1)`](./vsearch-udb2fasta.1.md),
[`vsearch-udbinfo(1)`](./vsearch-udbinfo.1.md),
[`vsearch-udb(5)`](../formats/vsearch-udb.5.md)


#(./fragments/footer.md)
