% vsearch-sintax(1) version 2.30.4 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./fragments/date.md)

# NAME

vsearch \-\-sintax --- classify sequences using the SINTAX algorithm


# SYNOPSIS

| **vsearch** **\-\-sintax** _fastxfile_ **\-\-db** _filename_ **\-\-tabbedout** _filename_ \[**\-\-sintax_cutoff** _real_] \[_options_]


# DESCRIPTION

The vsearch command `--sintax` classifies query sequences from a fasta or
fastq file using the SINTAX algorithm (Edgar 2016, doi:10.1101/074161), a
non-Bayesian method for taxonomic classification. The reference database is
specified with `--db`. Results are written to `--tabbedout`.

Classification works by finding the database sequence with the most shared
k-mers with each query across 100 bootstrap resamples of the query's k-mers.
For each bootstrap replicate, the taxonomy of the best-matching database
sequence is recorded. The frequency of agreement at each taxonomic rank across
the 100 replicates is reported as a bootstrap confidence value. Use
`--sintax_cutoff` to filter ranks below a confidence threshold.

The reference database must contain taxonomic annotations in the sequence
headers. Each header must include a `;tax=` field followed by a
comma-separated list of taxonomic identifiers. Each identifier starts with a
rank letter (`d` domain, `k` kingdom, `p` phylum, `c` class, `o` order,
`f` family, `g` genus, `s` species, `t` strain), a colon, and the taxon name.
Commas and semicolons are not allowed in taxon names. Example:

```text
>X80725;tax=d:Bacteria,p:Proteobacteria,c:Gammaproteobacteria,o:Enterobacteriales,f:Enterobacteriaceae,g:Escherichia,s:Escherichia_coli
ACGT...
```

The `--notrunclabels` option is **on by default** for `--sintax`, allowing
spaces in taxonomic identifiers.

When ties occur between database sequences with equally many k-mer matches,
the shortest (and then earliest) sequence is chosen by default. The option
`--sintax_random` is strongly recommended instead, as it breaks ties by a
random draw and avoids a bias towards shorter reference sequences.

For reproducible results with a fixed random seed, use `--randseed` together
with `--threads 1`. With multiple threads, sequences may be processed in
varying order across runs, making results non-reproducible even with a fixed
seed.

Both strands can be searched with `--strand both`. Databases in UDB format are
supported (see [`vsearch-udb(5)`](../formats/vsearch-udb.5.md)). This command
is multi-threaded.


# OPTIONS

## mandatory options

#(./fragments/option_db_sintax.md)

#(./fragments/option_tabbedout_sintax.md)


## core options

#(./fragments/option_dbmask.md)

#(./fragments/option_randseed.md)

#(./fragments/option_sintax_cutoff.md)

#(./fragments/option_sintax_random.md)

#(./fragments/option_strand.md)

#(./fragments/option_threads.md)


## secondary options

#(./fragments/option_bzip2_decompress.md)

#(./fragments/option_fastq_ascii.md)

#(./fragments/option_fastq_qmax.md)

#(./fragments/option_fastq_qmin.md)

#(./fragments/option_gzip_decompress.md)

#(./fragments/option_label_suffix.md)

#(./fragments/option_log.md)

#(./fragments/option_maxseqlength.md)

#(./fragments/option_minseqlength_32.md)

#(./fragments/option_no_progress.md)

#(./fragments/option_notrunclabels.md)

#(./fragments/option_quiet.md)

#(./fragments/option_wordlength.md)


# EXAMPLES

Classify query sequences against a taxonomy-annotated reference database,
writing all confidence values to the output:

```sh
vsearch \
    --sintax queries.fasta \
    --db reference.fasta \
    --tabbedout classification.tsv
```

Classify and report only ranks with at least 80% bootstrap support
(`--sintax_cutoff 0.8`), using the recommended tie-breaking option and
searching both strands:

```sh
vsearch \
    --sintax queries.fasta \
    --db reference.fasta \
    --sintax_cutoff 0.8 \
    --sintax_random \
    --strand both \
    --tabbedout classification.tsv
```

Produce reproducible results using a fixed random seed and a single thread:

```sh
vsearch \
    --sintax queries.fasta \
    --db reference.fasta \
    --sintax_cutoff 0.8 \
    --sintax_random \
    --randseed 42 \
    --threads 1 \
    --tabbedout classification.tsv
```


# SEE ALSO

[`vsearch-usearch_global(1)`](./vsearch-usearch_global.1.md),
[`vsearch-makeudb_usearch(1)`](./vsearch-makeudb_usearch.1.md),
[`vsearch-fasta(5)`](../formats/vsearch-fasta.5.md),
[`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md),
[`vsearch-udb(5)`](../formats/vsearch-udb.5.md)


#(./fragments/footer.md)
