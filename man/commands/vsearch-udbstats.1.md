% vsearch-udbstats(1) version 2.30.4 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./fragments/date.md)

# NAME

vsearch \-\-udbstats --- report statistics about indexed words in a UDB database file


# SYNOPSIS

| **vsearch** **\-\-udbstats** _udbfile_ \[_options_]


# DESCRIPTION

The vsearch command `--udbstats` reads a UDB database file and writes
statistics about its indexed words (*k*-mers) to the *standard output*
`stdout(3)`. The report gives the distribution of word frequencies
across the indexed sequences in the database.

See [`vsearch-udb(5)`](../formats/vsearch-udb.5.md) for a description
of the UDB file format.


# OPTIONS

## mandatory options

`--udbstats` *udbfile*
: Read the UDB database *udbfile* and report word statistics.


## secondary options

#(./fragments/option_log.md)

#(./fragments/option_no_progress.md)

#(./fragments/option_quiet.md)


## ignored options

#(./fragments/option_threads_not_multithreaded.md)


# EXAMPLES

Report word statistics for a UDB database:

```sh
vsearch \
    --udbstats db.udb
```


# SEE ALSO

[`vsearch-makeudb_usearch(1)`](./vsearch-makeudb_usearch.1.md),
[`vsearch-udb2fasta(1)`](./vsearch-udb2fasta.1.md),
[`vsearch-udbinfo(1)`](./vsearch-udbinfo.1.md),
[`vsearch-udb(5)`](../formats/vsearch-udb.5.md)


#(./fragments/footer.md)
