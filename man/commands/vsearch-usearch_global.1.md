% vsearch-usearch_global(1) version 2.30.4 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./fragments/date.md)

# NAME

vsearch \-\-usearch_global --- search sequences against a reference database


# SYNOPSIS

| **vsearch** **\-\-usearch_global** _fastxfile_ **\-\-db** _filename_ **\-\-id** _real_ (**\-\-alnout** | **\-\-biomout** | **\-\-blast6out** | **\-\-fastapairs** | **\-\-matched** | **\-\-mothur_shared_out** | **\-\-notmatched** | **\-\-otutabout** | **\-\-qsegout** | **\-\-samout** | **\-\-tsegout** | **\-\-uc** | **\-\-userout**) _filename_ \[_options_]


# DESCRIPTION

The vsearch command `--usearch_global` searches the query sequences in a
fasta or fastq file against a reference database (`--db`), using global
pairwise alignment (Needleman-Wunsch). The database can be a fasta or
fastq file, or a preformatted UDB database (see
[`vsearch-makeudb_usearch(1)`](./vsearch-makeudb_usearch.1.md)).

For each query, vsearch first pre-filters the database by counting shared
k-mers (words), then performs global pairwise alignments on the most
promising candidates. By default, the search stops after `--maxaccepts`
hits are accepted or `--maxrejects` candidates fail the identity threshold.
Setting both to 0 searches the entire database. Using values of `--id`
below 0.5 is unlikely to capture additional hits due to the k-mer
pre-filter.

The identity threshold is set with `--id`. By default, only the *plus*
strand of the query is compared to the database; use `--strand both` to
also check the reverse complement. Masking is applied with `--qmask`
(queries) and `--dbmask` (database).

At least one output option must be specified. This command is
multi-threaded.

To illustrate a search at 97% identity:

```text
Query file:    Database:       Results (--blast6out):

>q1            >t1             q1  t1  97.5  ...
ACGTACGT  -->  ACGTAGGT  -->   q2  *   *     ... (no hit)
>q2            >t2
TTTTTTTT       ACGTACGT
```


# OPTIONS

## mandatory options

`--db` and `--id` must both be specified, along with at least one output
option.

#(./fragments/option_db_usearch_global.md)

#(./fragments/option_id_search.md)


## core options

#(./fragments/option_dbmask.md)

#(./fragments/option_iddef.md)

#(./fragments/option_maxaccepts.md)

#(./fragments/option_maxrejects.md)

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

#(./fragments/option_idprefix.md)

#(./fragments/option_idsuffix.md)

#(./fragments/option_label_suffix.md)

#(./fragments/option_lca_cutoff.md)

#(./fragments/option_lcaout.md)

#(./fragments/option_leftjust.md)

#(./fragments/option_lengthout.md)

#(./fragments/option_log.md)

#(./fragments/option_matched.md)

#(./fragments/option_maxdiffs.md)

#(./fragments/option_maxgaps.md)

#(./fragments/option_maxhits.md)

#(./fragments/option_maxid.md)

#(./fragments/option_maxqsize.md)

#(./fragments/option_maxqt.md)

#(./fragments/option_maxseqlength.md)

#(./fragments/option_maxsizeratio.md)

#(./fragments/option_maxsl.md)

#(./fragments/option_maxsubs.md)

#(./fragments/option_mid.md)

#(./fragments/option_mincols.md)

#(./fragments/option_minqt.md)

#(./fragments/option_minseqlength_1.md)

#(./fragments/option_minsizeratio.md)

#(./fragments/option_minsl.md)

#(./fragments/option_mintsize.md)

#(./fragments/option_minwordmatches.md)

#(./fragments/option_mothur_shared_out.md)

#(./fragments/option_n_mismatch.md)

#(./fragments/option_no_progress.md)

#(./fragments/option_notmatched.md)

#(./fragments/option_notrunclabels.md)

#(./fragments/option_otutabout.md)

#(./fragments/option_output_no_hits.md)

#(./fragments/option_qsegout.md)

#(./fragments/option_query_cov.md)

#(./fragments/option_quiet.md)

#(./fragments/option_relabel.md)

#(./fragments/option_relabel_keep.md)

#(./fragments/option_relabel_md5.md)

#(./fragments/option_relabel_self.md)

#(./fragments/option_relabel_sha1.md)

#(./fragments/option_rightjust.md)

#(./fragments/option_rowlen.md)

#(./fragments/option_samheader.md)

#(./fragments/option_samout.md)

#(./fragments/option_sample.md)

#(./fragments/option_self.md)

#(./fragments/option_selfid.md)

#(./fragments/option_sizein.md)

#(./fragments/option_sizeout.md)

#(./fragments/option_target_cov.md)

#(./fragments/option_top_hits_only.md)

#(./fragments/option_tsegout.md)

#(./fragments/option_uc_search.md)

#(./fragments/option_uc_allhits.md)

#(./fragments/option_userfields.md)

#(./fragments/option_userout.md)

#(./fragments/option_weak_id.md)

#(./fragments/option_wordlength_8.md)

#(./fragments/option_xee.md)

#(./fragments/option_xlength.md)

#(./fragments/option_xsize.md)


## pairwise alignment options

These options modify the pairwise alignment scoring model. Modify with
caution.

#(./fragments/option_gapext.md)

#(./fragments/option_gapopen.md)

#(./fragments/option_match.md)

#(./fragments/option_mismatch.md)


## ignored options

These options are accepted for compatibility with usearch but have no
effect.

#(./fragments/option_band.md)

#(./fragments/option_fulldp.md)

#(./fragments/option_hspw.md)

#(./fragments/option_minhsp.md)

#(./fragments/option_pattern.md)

#(./fragments/option_slots.md)

#(./fragments/option_xdrop_nw.md)


# EXAMPLES

Search query sequences against a reference database at 97% identity and
write the results in BLAST-like tabular format:

```sh
vsearch \
    --usearch_global queries.fasta \
    --db reference.fasta \
    --id 0.97 \
    --blast6out results.tsv
```

Search both strands, report all hits per query (up to `--maxaccepts`), and
write matching query sequences to a fasta file:

```sh
vsearch \
    --usearch_global queries.fasta \
    --db reference.fasta \
    --id 0.97 \
    --strand both \
    --maxaccepts 5 \
    --matched matched.fasta \
    --notmatched unmatched.fasta
```

Classify query sequences against a taxonomically annotated database and
write last common ancestor assignments:

```sh
vsearch \
    --usearch_global queries.fasta \
    --db silva_tax.fasta \
    --id 0.97 \
    --maxaccepts 10 \
    --top_hits_only \
    --lcaout taxonomy.tsv \
    --alnout alignments.txt
```

Build an OTU table from multiple samples (sequences labelled with
`;sample=` annotations) against a set of OTU centroids:

```sh
vsearch \
    --usearch_global reads.fasta \
    --db otus.fasta \
    --id 0.97 \
    --otutabout otu_table.tsv
```


# SEE ALSO

[`vsearch-allpairs_global(1)`](./vsearch-allpairs_global.1.md),
[`vsearch-makeudb_usearch(1)`](./vsearch-makeudb_usearch.1.md),
[`vsearch-sintax(1)`](./vsearch-sintax.1.md),
[`vsearch-fasta(5)`](../formats/vsearch-fasta.5.md),
[`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md),
[`vsearch-udb(5)`](../formats/vsearch-udb.5.md),
[`vsearch-userfields(7)`](../misc/vsearch-userfields.7.md)


#(./fragments/footer.md)
