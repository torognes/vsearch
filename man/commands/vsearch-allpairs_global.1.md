% vsearch-allpairs_global(1) version 2.30.6 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./fragments/date.md)

# NAME

vsearch \-\-allpairs_global --- perform global pairwise alignments of all sequence pairs


# SYNOPSIS

| **vsearch** **\-\-allpairs_global** _fastafile_ (**\-\-acceptall** | **\-\-id** _real_) (**\-\-alnout** | **\-\-blast6out** | **\-\-fastapairs** | **\-\-matched** | **\-\-notmatched** | **\-\-qsegout** | **\-\-samout** | **\-\-tsegout** | **\-\-uc** | **\-\-userout**) _filename_ \[_options_]


# DESCRIPTION

The vsearch command `--allpairs_global` performs optimal global pairwise
alignments (Needleman-Wunsch) for all pairs of sequences in a fasta file.
Each sequence is compared to all sequences that follow it in the file,
for a total of n\*(n-1)/2 comparisons where n is the number of sequences.

Sequences are compared on their *plus* strand only. Either `--acceptall`
(to write all alignments to output) or `--id` (to set a minimum identity
threshold) must be specified. Most accept/reject options from the searching
section also apply.

Masking is applied as specified with `--qmask` and `--hardmask`.

At least one output option must be specified. This command is
multi-threaded.

To illustrate the comparisons performed on a three-sequence file:

```text
Input:    Pairs compared (n=3, total=3):

>s1       s1 vs s2
AAAA      s1 vs s3
>s2  -->  s2 vs s3
CCCC
>s3
TTTT
```


# OPTIONS

## mandatory options

Either `--acceptall` or `--id` must be specified, along with at least one
output option.

#(./fragments/option_acceptall.md)

#(./fragments/option_id_search.md)


## core options

#(./fragments/option_iddef.md)

#(./fragments/option_maxaccepts.md)

#(./fragments/option_maxrejects.md)

#(./fragments/option_qmask.md)

#(./fragments/option_threads.md)


## secondary options

#(./fragments/option_alnout.md)

#(./fragments/option_blast6out.md)

#(./fragments/option_bzip2_decompress.md)

#(./fragments/option_fasta_width.md)

#(./fragments/option_fastapairs.md)

#(./fragments/option_gzip_decompress.md)

#(./fragments/option_hardmask.md)

#(./fragments/option_idprefix.md)

#(./fragments/option_idsuffix.md)

#(./fragments/option_label_suffix.md)

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

#(./fragments/option_n_mismatch.md)

#(./fragments/option_no_progress.md)

#(./fragments/option_notmatched.md)

#(./fragments/option_notrunclabels.md)

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

#(./fragments/option_uc_allpairs.md)

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

Compare all pairs in a fasta file and write alignments above 97% identity
to a BLAST-like tabular file:

```sh
vsearch \
    --allpairs_global sequences.fasta \
    --id 0.97 \
    --blast6out results.tsv
```

Write all pairwise alignments regardless of identity, using a uclust-like
tabular format:

```sh
vsearch \
    --allpairs_global sequences.fasta \
    --acceptall \
    --uc results.uc
```

Compare all pairs and write matching query sequences to a fasta file,
using 8 threads:

```sh
vsearch \
    --allpairs_global sequences.fasta \
    --id 0.90 \
    --threads 8 \
    --matched matched.fasta \
    --notmatched unmatched.fasta
```


# SEE ALSO

[`vsearch-usearch_global(1)`](./vsearch-usearch_global.1.md),
[`vsearch-cigar(5)`](../formats/vsearch-cigar.5.md),
[`vsearch-fasta(5)`](../formats/vsearch-fasta.5.md),
[`vsearch-userfields(7)`](../misc/vsearch-userfields.7.md)


#(./fragments/footer.md)
