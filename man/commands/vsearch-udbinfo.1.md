% vsearch-udbinfo(1) version 2.30.6 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./fragments/date.md)

# NAME

vsearch \-\-udbinfo --- show information about a UDB database file


# SYNOPSIS

| **vsearch** **\-\-udbinfo** _udbfile_ \[_options_]


# DESCRIPTION

The vsearch command `--udbinfo` reads a UDB database file and writes
information about it to the *standard output* `stdout(3)`. The
information includes the number of sequences, the word length used for
the index, and masking settings.

See [`vsearch-udb(5)`](../formats/vsearch-udb.5.md) for a description
of the UDB file format.


# OPTIONS

## mandatory options

`--udbinfo` *udbfile*
: Read and inspect the UDB database *udbfile*. This option is
  mandatory.


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
