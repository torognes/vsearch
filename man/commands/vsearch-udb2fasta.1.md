% vsearch-udb2fasta(1) version 2.30.4 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./fragments/date.md)

# NAME

vsearch \-\-udb2fasta --- extract sequences from a UDB database file into a fasta file


# SYNOPSIS

| **vsearch** **\-\-udb2fasta** _udbfile_ **\-\-output** _fastafile_ \[_options_]


# DESCRIPTION

The vsearch command `--udb2fasta` reads a UDB database file and writes
its sequences to a fasta file. This is the reverse of
`--makeudb_usearch`.

Both `--udb2fasta` and `--output` must be specified.

See [`vsearch-udb(5)`](../formats/vsearch-udb.5.md) for a description
of the UDB file format.


# OPTIONS

## mandatory options

`--udb2fasta` *udbfile*
: Read sequences from the UDB database *udbfile*. This option is
  mandatory.

#(./fragments/option_output_udb2fasta.md)


## secondary options

#(./fragments/option_fasta_width.md)

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

#(./fragments/option_xee.md)

#(./fragments/option_xlength.md)

#(./fragments/option_xsize.md)


## ignored options

#(./fragments/option_threads_not_multithreaded.md)


# EXAMPLES

Extract sequences from a UDB file to fasta:

```sh
vsearch \
    --udb2fasta db.udb \
    --output db.fasta
```

Extract sequences and relabel them with a prefix:

```sh
vsearch \
    --udb2fasta db.udb \
    --relabel "seq" \
    --output relabelled.fasta
```


# SEE ALSO

[`vsearch-makeudb_usearch(1)`](./vsearch-makeudb_usearch.1.md),
[`vsearch-udbinfo(1)`](./vsearch-udbinfo.1.md),
[`vsearch-udbstats(1)`](./vsearch-udbstats.1.md),
[`vsearch-fasta(5)`](../formats/vsearch-fasta.5.md),
[`vsearch-udb(5)`](../formats/vsearch-udb.5.md)


#(./fragments/footer.md)
