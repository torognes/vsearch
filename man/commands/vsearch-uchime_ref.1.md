% vsearch-uchime_ref(1) version 2.31.0 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./fragments/date.md)

# NAME

vsearch \-\-uchime_ref --- detect chimeras using a reference database


# SYNOPSIS

| **vsearch** **\-\-uchime_ref** _fastafile_ (**\-\-chimeras** | **\-\-nonchimeras** | **\-\-uchimealns** | **\-\-uchimeout**) _filename_ **\-\-db** _dbfile_ \[_options_]


# DESCRIPTION

The vsearch command `--uchime_ref` detects chimeric sequences present
in the fasta-formatted *fastafile* by comparing them against a
reference database of chimera-free sequences (option `--db`).
Sequences are compared on their *plus* strand only; `--strand both`
is not supported by `--uchime_ref` and is rejected.

Chimera detection is based on a scoring function controlled by five
options: `--dn`, `--mindiffs`, `--mindiv`, `--minh`, and `--xn`. The
algorithm identifies candidate chimeras by finding three-way
alignments where a query sequence can be modelled as a mosaic of two
parent sequences from the reference database.

Chimeras can only be detected if their parents, or sufficiently close
relatives, are present in the reference database. Unlike the *de
novo* methods, `--uchime_ref` does not require abundance annotations.
Multithreading is supported.

Both `--db` and at least one output option must be specified.

See also `--uchime_denovo`, `--uchime2_denovo`, and
`--uchime3_denovo` for *de novo* chimera detection without a
reference database.


# OPTIONS

## mandatory options

`--uchime_ref` *fastafile*
: Detect chimeras in the fasta-formatted *fastafile* using a
  reference database. This option is mandatory.

#(./fragments/option_db_uchime_ref.md)

At least one of the following output options must be specified:

#(./fragments/option_chimeras.md)

#(./fragments/option_nonchimeras.md)

#(./fragments/option_uchimealns.md)

#(./fragments/option_uchimeout.md)

The `--borderline` option can also produce output, but only when
combined with at least one of the output options listed above:

#(./fragments/option_borderline.md)


## core options

#(./fragments/option_abskew_2.md)

#(./fragments/option_dbmask.md)

#(./fragments/option_dn.md)

#(./fragments/option_mindiffs.md)

#(./fragments/option_mindiv.md)

#(./fragments/option_minh.md)

#(./fragments/option_self_uchime_ref.md)

#(./fragments/option_selfid_uchime_ref.md)

#(./fragments/option_threads.md)

#(./fragments/option_uchimeout5.md)

#(./fragments/option_xn.md)


## secondary options

#(./fragments/option_alignwidth_80.md)

#(./fragments/option_fasta_score.md)

#(./fragments/option_fasta_width.md)

#(./fragments/option_hardmask.md)

#(./fragments/option_label_suffix.md)

#(./fragments/option_lengthout.md)

#(./fragments/option_log.md)

#(./fragments/option_maxseqlength.md)

#(./fragments/option_minseqlength_1.md)

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

#(./fragments/option_sizeout.md)

#(./fragments/option_xee.md)

#(./fragments/option_xlength.md)

#(./fragments/option_xsize.md)


## pairwise alignment options

These options modify the parameters of the pairwise alignment model.
Modify with caution.

#(./fragments/option_gapext.md)

#(./fragments/option_gapopen.md)

#(./fragments/option_match.md)

#(./fragments/option_mismatch.md)


## ignored options

#(./fragments/option_sizein.md)
: Ignored by `--uchime_ref`.


# EXAMPLES

Detect chimeras using a reference database and write non-chimeric
sequences to a file:

```sh
vsearch \
    --uchime_ref amplicons.fasta \
    --db silva_138_db.fasta \
    --nonchimeras clean.fasta
```

Use multiple threads and write a detailed tab-separated report:

```sh
vsearch \
    --uchime_ref amplicons.fasta \
    --db silva_138_db.fasta \
    --threads 4 \
    --nonchimeras clean.fasta \
    --uchimeout chimera_report.tsv
```

Write all three output categories (chimeras, non-chimeras and
borderline entries):

```sh
vsearch \
    --uchime_ref amplicons.fasta \
    --db silva_138_db.fasta \
    --nonchimeras clean.fasta \
    --chimeras chimeras.fasta \
    --borderline borderline.fasta
```


# SEE ALSO

[`vsearch-uchime_denovo(1)`](./vsearch-uchime_denovo.1.md),
[`vsearch-uchime2_denovo(1)`](./vsearch-uchime2_denovo.1.md),
[`vsearch-uchime3_denovo(1)`](./vsearch-uchime3_denovo.1.md),
[`vsearch-chimeras_denovo(1)`](./vsearch-chimeras_denovo.1.md),
[`vsearch-fasta(5)`](../formats/vsearch-fasta.5.md),
[`vsearch-udb(5)`](../formats/vsearch-udb.5.md)


#(./fragments/footer.md)
