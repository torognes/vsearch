% vsearch-cluster_smallmem(1) version 2.30.4 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./fragments/date.md)

# NAME

vsearch \-\-cluster_smallmem --- clusterize pre-sorted sequences


# SYNOPSIS

| **vsearch** **\-\-cluster_smallmem** _fastafile_ **\-\-id** _real_ (**\-\-alnout** | **\-\-biomout** | **\-\-blast6out** | **\-\-centroids** | **\-\-clusters** | **\-\-mothur_shared_out** | **\-\-msaout** | **\-\-otutabout** | **\-\-profile** | **\-\-samout** | **\-\-uc** | **\-\-userout**) _filename_ \[_options_]


# DESCRIPTION

The vsearch command `--cluster_smallmem` groups the fasta sequences in
*fastafile* into clusters using a greedy, heuristic, centroid-based
algorithm. Unlike `--cluster_fast` and `--cluster_size`, this command
does not sort the input sequences before clustering; sequences are
expected to be sorted by decreasing length. Use `--usersort` to
override this requirement and allow sequences in any order. At least
one output option must be specified.

For each query sequence (in input order), vsearch compares it to all
existing cluster centroids. If the query is similar enough to a
centroid (as determined by `--id`), it is assigned to that cluster;
otherwise, a new cluster is seeded with the query as its centroid.

Sequences are compared using global pairwise alignment
(Needleman-Wunsch). The number of comparisons is limited by
`--maxaccepts` and `--maxrejects`.

`--cluster_fast` (see
[`vsearch-cluster_fast(1)`](./vsearch-cluster_fast.1.md)) and
`--cluster_size` (see
[`vsearch-cluster_size(1)`](./vsearch-cluster_size.1.md)) perform the
same clustering but automatically sort the input sequences beforehand.

See [`vsearch-fasta(5)`](../formats/vsearch-fasta.5.md) for a
description of the input format.


# OPTIONS

## mandatory options

`--cluster_smallmem` *fastafile*
: Read and clusterize fasta sequences from *fastafile*. Sequences must
  be sorted by decreasing length, unless `--usersort` is used. This
  option is mandatory.

#(./fragments/option_id.md)


## core options

#(./fragments/option_centroids.md)

#(./fragments/option_clusterout_id.md)

#(./fragments/option_clusterout_sort.md)

#(./fragments/option_clusters.md)

#(./fragments/option_consout.md)

#(./fragments/option_iddef.md)

#(./fragments/option_maxaccepts.md)

#(./fragments/option_maxrejects.md)

#(./fragments/option_msaout.md)

#(./fragments/option_profile.md)

#(./fragments/option_qmask.md)

#(./fragments/option_sizein.md)

#(./fragments/option_sizeorder.md)

#(./fragments/option_sizeout.md)

#(./fragments/option_strand.md)

#(./fragments/option_uc.md)

#(./fragments/option_usersort.md)


## secondary options

#(./fragments/option_alnout.md)

#(./fragments/option_biomout.md)

#(./fragments/option_blast6out.md)

#(./fragments/option_bzip2_decompress.md)

#(./fragments/option_fastapairs.md)

#(./fragments/option_fasta_width.md)

#(./fragments/option_gapext.md)

#(./fragments/option_gapopen.md)

#(./fragments/option_gzip_decompress.md)

#(./fragments/option_hardmask.md)

#(./fragments/option_idprefix.md)

#(./fragments/option_idsuffix.md)

#(./fragments/option_label_suffix.md)

#(./fragments/option_leftjust.md)

#(./fragments/option_lengthout.md)

#(./fragments/option_log.md)

#(./fragments/option_match.md)

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

#(./fragments/option_minseqlength_32.md)

#(./fragments/option_minqt.md)

#(./fragments/option_minsizeratio.md)

#(./fragments/option_minsl.md)

#(./fragments/option_mintsize.md)

#(./fragments/option_minwordmatches.md)

#(./fragments/option_mismatch.md)

#(./fragments/option_mothur_shared_out.md)

#(./fragments/option_n_mismatch.md)

#(./fragments/option_no_progress.md)

#(./fragments/option_notmatched.md)

#(./fragments/option_notrunclabels.md)

#(./fragments/option_otutabout.md)

#(./fragments/option_output_no_hits.md)

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

#(./fragments/option_qsegout.md)

#(./fragments/option_target_cov.md)

#(./fragments/option_threads.md)

#(./fragments/option_top_hits_only.md)

#(./fragments/option_tsegout.md)

#(./fragments/option_userfields.md)

#(./fragments/option_userout.md)

#(./fragments/option_weak_id.md)

#(./fragments/option_wordlength.md)

#(./fragments/option_xlength.md)

#(./fragments/option_xsize.md)


## ignored options

#(./fragments/option_band.md)

#(./fragments/option_cons_truncate.md)

#(./fragments/option_fulldp.md)

#(./fragments/option_hspw.md)

#(./fragments/option_minhsp.md)

#(./fragments/option_pattern.md)

#(./fragments/option_slots.md)

#(./fragments/option_xdrop_nw.md)


# EXAMPLES

Cluster sequences at 97% identity, assuming the input is already
sorted by decreasing length:

```sh
vsearch \
    --cluster_smallmem sequences.fasta \
    --id 0.97 \
    --sizein \
    --sizeout \
    --centroids centroids.fasta
```

Cluster sequences in any order by using `--usersort`:

```sh
vsearch \
    --cluster_smallmem sequences.fasta \
    --id 0.97 \
    --usersort \
    --sizein \
    --sizeout \
    --centroids centroids.fasta
```


# SEE ALSO

[`vsearch-cluster_fast(1)`](./vsearch-cluster_fast.1.md),
[`vsearch-cluster_size(1)`](./vsearch-cluster_size.1.md),
[`vsearch-cluster_unoise(1)`](./vsearch-cluster_unoise.1.md),
[`vsearch-derep_fulllength(1)`](./vsearch-derep_fulllength.1.md),
[`vsearch-fastx_uniques(1)`](./vsearch-fastx_uniques.1.md),
[`vsearch-uchime_denovo(1)`](./vsearch-uchime_denovo.1.md),
[`vsearch-fasta(5)`](../formats/vsearch-fasta.5.md)


#(./fragments/footer.md)
