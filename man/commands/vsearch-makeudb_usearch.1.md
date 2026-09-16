% vsearch-makeudb_usearch(1) version 2.32.0 | vsearch manual
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

The UDB file holds one 32-bit counter per possible *k*-mer, whether or
not the database contains that *k*-mer. Those 4^*wordlength* counters
are 256 kB at the default word length of 8, but 1 GB at word length 14
and 4 GB at 15, whatever the size of *fastxfile*, so raising
`--wordlength` grows the output file as well as the memory needed to
build it (see `--wordlength`).

Of the work this command performs, only DUST masking is distributed
over several threads. Reading the input, building the *k*-mer index and
writing the UDB file are single-threaded, and so is the letter
replacement `--hardmask` asks for. `--threads` therefore shortens the
masking step and nothing else: it helps with the default
`--dbmask dust`, and has no measurable effect with `--dbmask none` or
`--dbmask soft`. What it can save is bounded by the share of the run
spent masking, which grows with the length of the sequences, and on
short ones the gain stops increasing beyond a handful of threads.

See [`vsearch-udb(5)`](../formats/vsearch-udb.5.md) for a description
of the UDB file format.


# OPTIONS

## mandatory options

`--makeudb_usearch` *fastxfile*
: Read fasta or fastq sequences from *fastxfile* and create a UDB
  database (quality values are ignored).

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

#(./fragments/option_threads.md)


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
