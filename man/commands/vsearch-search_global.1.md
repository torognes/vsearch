% vsearch-search_global(1) version 2.33.0 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./fragments/date.md)

# NAME

vsearch \-\-search_global --- exhaustively search sequences against a reference database


# SYNOPSIS

| **vsearch** **\-\-search_global** _fastxfile_ **\-\-db** _filename_ **\-\-id** _real_ (**\-\-alnout** | **\-\-biomout** | **\-\-blast6out** | **\-\-dbmatched** | **\-\-dbnotmatched** | **\-\-fastapairs** | **\-\-lcaout** | **\-\-matched** | **\-\-mothur_shared_out** | **\-\-notmatched** | **\-\-otutabout** | **\-\-qsegout** | **\-\-samout** | **\-\-tsegout** | **\-\-uc** | **\-\-userout**) _filename_ \[_options_]


# DESCRIPTION

The vsearch command `--search_global` searches the query sequences in a
fasta or fastq file against a reference database (`--db`), using global
pairwise alignment (Needleman-Wunsch). It is the exhaustive counterpart
of [`vsearch-usearch_global(1)`](./vsearch-usearch_global.1.md): **every
database sequence is aligned against every query**. No *k*-mer index is
built, no pre-filter selects candidates, and the search does not stop
early. A hit is therefore found whatever its identity, and the results do
not depend on the order the database happens to be in.

The database can be a fasta or fastq file, or a preformatted UDB database
(see [`vsearch-makeudb_usearch(1)`](./vsearch-makeudb_usearch.1.md)); of a
UDB, only the sequences are read, the index it carries being of no use
here.

This is slow. Expect a runtime proportional to the number of queries
times the size of the database: as an order of magnitude, about five
seconds per query against 200,000 references of 300 nucleotides on one
core. The command is multi-threaded, and the work divides cleanly, but no
thread count turns an exhaustive search into a quick one. Use
`--usearch_global` for routine work, and `--search_global` when you need
the guarantee that nothing was skipped --- typically on a small database,
or at an identity threshold too low for the pre-filter to respect.

**Mind the volume of output.** `--maxhits` defaults to 0, meaning
unlimited, and a low `--id` accepts nearly everything, so a run can write
one line *per database sequence and per query*: a thousand queries against
a database of 200,000 references is 2 x 10^8 lines. Set `--maxhits`, or
`--top_hits_only`, unless you really want them all.

**Mind what a low identity means.** Global alignment of two unrelated
sequences of similar length routinely reports 30% to 45% identity, because
the aligner still has to line them up end to end. Such hits are the
command working as designed, not a defect, and they are exactly the hits
`--usearch_global` never shows you. An identity threshold below roughly
0.5 will return them in quantity (see also
[`vsearch-pairwise_alignment_parameters(7)`](../misc/vsearch-pairwise_alignment_parameters.7.md)).

Alignment is global, not local: a query is aligned end-to-end against a
target rather than searched for as a subsequence, and at most one
alignment is reported per database sequence and per strand. A motif
occurring several times within one long target therefore yields a single
hit, and two at most with `--strand both`. Enumerating every occurrence
calls for a local-alignment tool instead.

The identity threshold is set with `--id`. By default, only the *plus*
strand of the query is compared to the database; use `--strand both` to
also check the reverse complement. Masking is applied with `--qmask`
(queries) and `--dbmask` (database); with soft masking it has no effect on
the alignment here, since masking only ever steered the *k*-mer index, but
together with `--hardmask` it replaces the masked bases with Ns and does
change the result.

At least one output option must be specified. This command is
multi-threaded: the queries are distributed over the available threads,
so the order of the entries written to `--alnout`, `--blast6out`,
`--fastapairs`, `--lcaout`, `--matched`, `--notmatched`, `--qsegout`,
`--samout`, `--tsegout`, `--uc` and `--userout` may vary from run to run
when more than one thread is used. The `--biomout`, `--dbmatched`,
`--dbnotmatched`, `--mothur_shared_out` and `--otutabout` tables are
assembled after the search, or written in database order, and keep a
stable order. The results themselves do not depend on the thread count.

To illustrate a search at 97% identity, where the second query has no
target above the threshold but is still compared to every one of them:

```text
Query file:    Database:       Results (--blast6out):

>q1            >t1             q1  t1  97.5  ...
ACGTACGT  -->  ACGTAGGT  -->   q2  *   *     ... (no hit; this line is
>q2            >t2                 only written with --output_no_hits)
TTTTTTTT       ACGTACGT
```


# OPTIONS

## mandatory options

`--db` and `--id` must both be specified, along with at least one output
option.

#(./fragments/option_db_usearch_global.md)

#(./fragments/option_id_search_global.md)


## core options

#(./fragments/option_dbmask.md)

#(./fragments/option_iddef.md)

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

#(./fragments/option_xee.md)

#(./fragments/option_xlength.md)

#(./fragments/option_xsize.md)


## pairwise alignment options

#(./fragments/option_gapext.md)

#(./fragments/option_gapopen.md)

#(./fragments/option_match.md)

#(./fragments/option_mismatch.md)


## ignored options

These options are accepted for compatibility with usearch, or so that a
command line written for `--usearch_global` can be reused unchanged, but
have no effect.

#(./fragments/option_band.md)

#(./fragments/option_fulldp.md)

#(./fragments/option_hspw.md)

#(./fragments/option_maxaccepts_ignored_search_global.md)

#(./fragments/option_maxrejects_ignored_search_global.md)

#(./fragments/option_minhsp.md)

#(./fragments/option_minwordmatches_ignored_search_global.md)

#(./fragments/option_pattern.md)

#(./fragments/option_slots.md)

#(./fragments/option_wordlength_ignored_search_global.md)

#(./fragments/option_xdrop_nw.md)


# EXAMPLES

Search a small set of queries against a reference database, reporting
every hit above 50% identity:

```sh
vsearch \
    --search_global queries.fasta \
    --db references.fasta \
    --id 0.5 \
    --blast6out hits.tsv
```

Keep only the best hit per query, which is the usual way to tame the
output of an exhaustive search:

```sh
vsearch \
    --search_global queries.fasta \
    --db references.fasta \
    --id 0.5 \
    --maxhits 1 \
    --blast6out best_hit.tsv
```

Find the remote relatives of a handful of sequences, at an identity too
low for the word pre-filter of `--usearch_global` to respect, and write
the alignments so they can be inspected:

```sh
vsearch \
    --search_global orphans.fasta \
    --db references.fasta \
    --id 0.4 \
    --maxhits 10 \
    --alnout alignments.txt
```

Search both strands, against a preformatted UDB database (only its
sequences are read):

```sh
vsearch \
    --search_global queries.fasta \
    --db references.udb \
    --id 0.7 \
    --strand both \
    --userout hits.tsv \
    --userfields query+target+id+qstrand
```


# SEE ALSO

[`vsearch-usearch_global(1)`](./vsearch-usearch_global.1.md),
[`vsearch-allpairs_global(1)`](./vsearch-allpairs_global.1.md),
[`vsearch-search_exact(1)`](./vsearch-search_exact.1.md),
[`vsearch-makeudb_usearch(1)`](./vsearch-makeudb_usearch.1.md),
[`vsearch-cigar(5)`](../formats/vsearch-cigar.5.md),
[`vsearch-fasta(5)`](../formats/vsearch-fasta.5.md),
[`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md),
[`vsearch-udb(5)`](../formats/vsearch-udb.5.md),
[`vsearch-pairwise_alignment_parameters(7)`](../misc/vsearch-pairwise_alignment_parameters.7.md),
[`vsearch-usearch(7)`](../misc/vsearch-usearch.7.md),
[`vsearch-userfields(7)`](../misc/vsearch-userfields.7.md)


#(./fragments/footer.md)
