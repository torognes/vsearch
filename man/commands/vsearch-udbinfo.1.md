% vsearch-udbinfo(1) version 2.31.0 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./fragments/date.md)

# NAME

vsearch \-\-udbinfo --- show information about a UDB database file


# SYNOPSIS

| **vsearch** **\-\-udbinfo** _udbfile_ \[_options_]


# DESCRIPTION

The vsearch command `--udbinfo` reads a UDB database file and writes
information about it to the *standard error* `stderr(3)`. The
information includes the number of sequences (Seqs), the number of
bits per sequence index entry (SeqIx bits), the alphabet (Alpha), the
word length used for the index (Word width), the number of index slots
(Slots), the dictionary size (Dict size), and the DBstep and DBAccel
settings. Masking information is not stored in UDB files and is
therefore not reported.

Only Seqs, Word width, DBstep and DBAccel are read from the file;
the rest are fixed by the format. In particular Slots is header field
`buffer[11]`, which `--makeudb_usearch` writes as zero (see
[`vsearch-udb(5)`](../formats/vsearch-udb.5.md)), so it reads `0` for
every UDB file vsearch has written. Dict size is the real number of
index slots, 4 raised to the power of the word length, and
[`vsearch-udbstats(1)`](./vsearch-udbstats.1.md) reports the same value
under the name Slots.

Because `--udbinfo` reads only the file's 200-byte header, it does not
detect a corrupt or truncated body; the other UDB commands do. Note
also that `--quiet` suppresses the report entirely, since `stderr` is
where it goes: combine `--quiet` with `--log` to send it to a file
instead.

See [`vsearch-udb(5)`](../formats/vsearch-udb.5.md) for a description
of the UDB file format.


# OPTIONS

## mandatory options

`--udbinfo` *udbfile*
: Read and inspect the UDB database *udbfile*. This option is
  mandatory. As UDB files cannot be read from pipes, *udbfile* must be a
  seekable file path (see
  [`vsearch-udb(5)`](../formats/vsearch-udb.5.md)); a pipe is rejected
  with an error.


## secondary options

#(./fragments/option_log.md)

#(./fragments/option_quiet.md)


## ignored options

#(./fragments/option_threads_not_multithreaded.md)


# EXAMPLES

Display information about a UDB database:

```sh
vsearch \
    --udbinfo db.udb
```


# SEE ALSO

[`vsearch-makeudb_usearch(1)`](./vsearch-makeudb_usearch.1.md),
[`vsearch-udb2fasta(1)`](./vsearch-udb2fasta.1.md),
[`vsearch-udbstats(1)`](./vsearch-udbstats.1.md),
[`vsearch-udb(5)`](../formats/vsearch-udb.5.md)


#(./fragments/footer.md)
