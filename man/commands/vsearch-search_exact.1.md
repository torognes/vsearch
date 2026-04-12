% vsearch-search_exact(1) version 2.30.4 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./fragments/date.md)

# NAME

vsearch \-\-search_exact --- search for exact full-length matches against a database


# SYNOPSIS

| **vsearch** **\-\-search_exact** _fastxfile_ **\-\-db** _filename_ (**\-\-alnout** | **\-\-biomout** | **\-\-blast6out** | **\-\-fastapairs** | **\-\-matched** | **\-\-mothur_shared_out** | **\-\-notmatched** | **\-\-otutabout** | **\-\-qsegout** | **\-\-samout** | **\-\-tsegout** | **\-\-uc** | **\-\-userout**) _filename_ \[_options_]


# DESCRIPTION

The vsearch command `--search_exact` searches the query sequences in a
fasta or fastq file against a database of target sequences (`--db`),
reporting only 100% exact full-length matches. It is much faster than
`--usearch_global` for this use case.

Because only exact matches are reported, `--id`, `--maxaccepts`, and
`--maxrejects` do not apply and are not accepted. Masking is controlled
with `--qmask` (queries) and `--dbmask` (database). By default only the
*plus* strand is searched; use `--strand both` to also check the reverse
complement.

At least one output option must be specified. This command is
multi-threaded.


# OPTIONS

## mandatory options

`--db` must be specified, along with at least one output option.

#(./fragments/option_db_search_exact.md)


## core options

#(./fragments/option_dbmask.md)

#(./fragments/option_qmask.md)

#(./fragments/option_strand.md)

#(./fragments/option_threads.md)


## secondary options

#(./fragments/option_alnout.md)

#(./fragments/option_biomout.md)

#(./fragments/option_blast6out.md)

#(./fragments/option_bzip2_decompress.md)

#(./fragments/option_dbmatched.md)

#(./fragments/option_dbnotmatched.md)

#(./fragments/option_fasta_width.md)

#(./fragments/option_fastapairs.md)

#(./fragments/option_gzip_decompress.md)

#(./fragments/option_hardmask.md)

#(./fragments/option_label_suffix.md)

#(./fragments/option_lca_cutoff.md)

#(./fragments/option_lcaout.md)

#(./fragments/option_lengthout.md)

#(./fragments/option_log.md)

#(./fragments/option_matched.md)

#(./fragments/option_maxhits.md)

#(./fragments/option_maxqsize.md)

#(./fragments/option_maxqt.md)

#(./fragments/option_maxseqlength.md)

#(./fragments/option_maxsizeratio.md)

#(./fragments/option_maxsl.md)

#(./fragments/option_mincols.md)

#(./fragments/option_minqt.md)

#(./fragments/option_minseqlength_1.md)

#(./fragments/option_minsizeratio.md)

#(./fragments/option_minsl.md)

#(./fragments/option_mintsize.md)

#(./fragments/option_mothur_shared_out.md)

#(./fragments/option_no_progress.md)

#(./fragments/option_notmatched.md)

#(./fragments/option_notrunclabels.md)

#(./fragments/option_otutabout.md)

#(./fragments/option_output_no_hits.md)

#(./fragments/option_qsegout.md)

#(./fragments/option_quiet.md)

#(./fragments/option_relabel.md)

#(./fragments/option_relabel_keep.md)

#(./fragments/option_relabel_md5.md)

#(./fragments/option_relabel_self.md)

#(./fragments/option_relabel_sha1.md)

#(./fragments/option_rowlen.md)

#(./fragments/option_samheader.md)

#(./fragments/option_samout.md)

#(./fragments/option_sample.md)

#(./fragments/option_self.md)

#(./fragments/option_sizein.md)

#(./fragments/option_sizeout.md)

#(./fragments/option_top_hits_only.md)

#(./fragments/option_tsegout.md)

#(./fragments/option_uc_search.md)

#(./fragments/option_uc_allhits.md)

#(./fragments/option_userfields.md)

#(./fragments/option_userout.md)

#(./fragments/option_xee.md)

#(./fragments/option_xlength.md)

#(./fragments/option_xsize.md)


## pairwise alignment options

These options modify the pairwise alignment scoring model used when
writing alignment output (e.g., `--alnout`). Modify with caution.

#(./fragments/option_match.md)

#(./fragments/option_mismatch.md)


# EXAMPLES

Find all exact matches in a database and write results in BLAST-like
tabular format:

```sh
vsearch \
    --search_exact queries.fasta \
    --db reference.fasta \
    --blast6out results.tsv
```

Search both strands and write matching and non-matching queries to
separate fasta files:

```sh
vsearch \
    --search_exact queries.fasta \
    --db reference.fasta \
    --strand both \
    --matched exact_hits.fasta \
    --notmatched no_hits.fasta
```

Build an OTU table by mapping dereplicated reads to OTU centroids using
exact matching:

```sh
vsearch \
    --search_exact reads.fasta \
    --db otus.fasta \
    --otutabout otu_table.tsv \
    --threads 8
```


# SEE ALSO

[`vsearch-usearch_global(1)`](./vsearch-usearch_global.1.md),
[`vsearch-fasta(5)`](../formats/vsearch-fasta.5.md),
[`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md)


#(./fragments/footer.md)
