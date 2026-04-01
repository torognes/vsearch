% vsearch-rereplicate(1) version 2.30.6 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./fragments/date.md)

# NAME

vsearch \-\-rereplicate --- use abundance values to rereplicate fasta sequences


# SYNOPSIS

| **vsearch** **\-\-rereplicate** *fastafile* \-\-output *outputfile* \[*options*]


# DESCRIPTION

The vsearch command `--rereplicate` uses abundance values contained in
fasta header annotations (`;size=n`) to rereplicate sequences. An
input fasta entry is rereplicated *n* times if its abudance is *n*
(`;size=n`) (option `--sizein` is always implied). If abundance
annotations are missing, a warning is emitted. The label of the output
sequence remains the same, unless `--relabel`, `--relabel_self`,
`--relabel_sha1` or `--relabel_md5` are used. Output is written to the
file specified with the `--output` option, in fasta format. The output
file does not contain abundance information, unless `--sizeout` is
specified, in which case an abundance of 1 is used (`;size=1`).

The inverse operation is named *dereplication* and is performed by the
command `--fastx_uniques` (see
[`vsearch-fastx_uniques(1)`](./vsearch-fastx_uniques.1.md)).


# OPTIONS

## mandatory options

#(./fragments/option_output.md)


## core options

`--sizein`
: Use the abundance annotations (`[>@;]size=integer[;]`) present in
  sequence headers when reading the input fasta file. Entries without
  abundance annotations are assumed to be of `;size=1`, and a warning
  is emitted. This option is always implied when using `--rereplicate`.

`--sizeout`
: Add abundance annotation `;size=1` to sequence headers when writing
  the output fasta file.


## secondary options

#(./fragments/option_bzip2_decompress.md)

#(./fragments/option_fasta_width.md)

#(./fragments/option_gzip_decompress.md)

#(./fragments/option_label_suffix.md)

#(./fragments/option_lengthout.md)

#(./fragments/option_log.md)

#(./fragments/option_no_progress.md)

#(./fragments/option_notrunclabels.md)

#(./fragments/option_quiet.md)

#(./fragments/option_relabel.md)

#(./fragments/option_relabel_keep.md)

#(./fragments/option_relabel_md5.md)

#(./fragments/option_relabel_self.md)

#(./fragments/option_relabel_sha1.md)

#(./fragments/option_sample.md)

#(./fragments/option_xee.md)

#(./fragments/option_xlength.md)

#(./fragments/option_xsize.md)


## ignored options

#(./fragments/option_threads_not_multithreaded.md)


# EXAMPLES

Rereplicate the sequences in *dereplicated.fasta*. Write the
duplicated sequences to *bigfile.fasta*, in fasta format:

```sh
vsearch \
    --rereplicate dereplicated.fasta \
    --output bigfile.fasta
```


# SEE ALSO

[`vsearch-fastx_uniques(1)`](./vsearch-fastx_uniques.1.md),
[`vsearch-fasta(5)`](./vsearch-fasta.5.md)


#(./fragments/footer.md)
