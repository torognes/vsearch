% vsearch-makeudb_usearch(1) version 2.31.0 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./fragments/date.md)

# NAME

vsearch \-\-makeudb_usearch --- create a UDB database file from a fasta
or fastq file


# SYNOPSIS

| **vsearch** **\-\-makeudb_usearch** _fastxfile_ **\-\-output** _dbfile_ \[_options_]


# DESCRIPTION

The vsearch command `--makeudb_usearch` creates a UDB database file
from the sequences in *fastxfile*, in fasta or fastq format (quality
values are ignored and are not stored in the database). The UDB file is a
binary format that contains the sequences together with a *k*-mer
index, and can be quickly loaded into memory. Using a UDB file avoids
re-indexing the database on every run, which is worthwhile when the
same database is searched repeatedly with `--usearch_global` or
`--sintax`.

Both `--makeudb_usearch` and `--output` must be specified.

The database must contain at least one sequence. vsearch reports a fatal
error and stops if *fastxfile* is empty, or if all of its sequences are
discarded by `--minseqlength` (32 nucleotides by default for this
command), since a UDB file recording no sequence cannot be read back.

See [`vsearch-udb(5)`](../formats/vsearch-udb.5.md) for a description
of the UDB file format.


# OPTIONS

## mandatory options

`--makeudb_usearch` *fastxfile*
: Read fasta or fastq sequences from *fastxfile* and create a UDB
  database (quality values are ignored). This option is mandatory.

#(./fragments/option_output_makeudb.md)


## core options

#(./fragments/option_dbmask.md)

#(./fragments/option_hardmask.md)

#(./fragments/option_wordlength_8.md)


## secondary options

#(./fragments/option_bzip2_decompress.md)

#(./fragments/option_gzip_decompress.md)

#(./fragments/option_log.md)

#(./fragments/option_maxseqlength.md)

#(./fragments/option_minseqlength_32.md)

#(./fragments/option_no_progress.md)

#(./fragments/option_notrunclabels.md)

#(./fragments/option_quiet.md)


## ignored options

#(./fragments/option_threads_not_multithreaded.md)


# EXAMPLES

Create a UDB database from a fasta file with default settings:

```sh
vsearch \
    --makeudb_usearch db.fasta \
    --output db.udb
```

Create a UDB database without masking:

```sh
vsearch \
    --makeudb_usearch db.fasta \
    --dbmask none \
    --output db.udb
```

Use the resulting UDB file with `--usearch_global`:

```sh
vsearch \
    --makeudb_usearch db.fasta \
    --output db.udb

vsearch \
    --usearch_global queries.fasta \
    --db db.udb \
    --id 0.97 \
    --blast6out results.b6
```


# SEE ALSO

[`vsearch-udb2fasta(1)`](./vsearch-udb2fasta.1.md),
[`vsearch-udbinfo(1)`](./vsearch-udbinfo.1.md),
[`vsearch-udbstats(1)`](./vsearch-udbstats.1.md),
[`vsearch-fasta(5)`](../formats/vsearch-fasta.5.md),
[`vsearch-udb(5)`](../formats/vsearch-udb.5.md)


#(./fragments/footer.md)
