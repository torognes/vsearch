% vsearch-cluster_unoise(1) version 2.31.0 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./fragments/date.md)

# NAME

vsearch \-\-cluster_unoise --- denoise amplicon sequences using the UNOISE3 algorithm


# SYNOPSIS

| **vsearch** **\-\-cluster_unoise** _fastafile_ (**\-\-alnout** | **\-\-biomout** | **\-\-blast6out** | **\-\-centroids** | **\-\-clusters** | **\-\-mothur_shared_out** | **\-\-msaout** | **\-\-otutabout** | **\-\-profile** | **\-\-samout** | **\-\-uc** | **\-\-userout**) _filename_ \[_options_]


# DESCRIPTION

The vsearch command `--cluster_unoise` performs denoising of fasta
amplicon sequences according to the UNOISE version 3 algorithm by
Robert Edgar, but without the *de novo* chimera removal step. The
chimera removal step may be performed separately with
`--uchime3_denovo` (see
[`vsearch-uchime3_denovo(1)`](./vsearch-uchime3_denovo.1.md)).

The input sequences are expected to be sorted by decreasing abundance,
and must carry abundance annotations in their headers (e.g.
`;size=integer;`). Use `--sizein` to read this information.

At least one output option must be specified.

Unlike the standard clustering commands, `--cluster_unoise` does not
use a fixed identity threshold (`--id`). Instead, a sequence is
assigned to an existing cluster if two conditions are met: (1) the
sequence identity with the centroid is high enough, and (2) the
abundance ratio (the abundance of the candidate divided by the
abundance of the centroid) does not exceed *beta*. Beta is computed as
2^(-1 - alpha × distance), where *alpha* is set with `--unoise_alpha`
(default 2.0) and *distance* is the number of mismatches in the
alignment, ignoring gaps. As a result, sequences must be exponentially
less abundant as their distance from the centroid increases in order to
be grouped together; more abundant sequences at any distance will form
their own new clusters.

The minimum abundance threshold for input sequences is controlled by
`--minsize` (default 8). Sequences with an abundance below this value
are excluded before denoising.


# OPTIONS

## mandatory options

`--cluster_unoise` *fastafile*
: Read and denoise fasta amplicon sequences from *fastafile*. This
  option is mandatory.


## core options

#(./fragments/option_centroids.md)

#(./fragments/option_clusterout_id.md)

#(./fragments/option_clusterout_sort.md)

#(./fragments/option_clusters.md)

#(./fragments/option_consout.md)

#(./fragments/option_iddef.md)

#(./fragments/option_maxaccepts.md)

#(./fragments/option_maxrejects.md)

#(./fragments/option_minsize.md)

#(./fragments/option_msaout.md)

#(./fragments/option_profile.md)

#(./fragments/option_qmask.md)

#(./fragments/option_sizein.md)

#(./fragments/option_sizeout.md)

#(./fragments/option_strand.md)

#(./fragments/option_uc.md)

#(./fragments/option_unoise_alpha.md)


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

#(./fragments/option_id.md)

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

Dereplicate sequences before denoising, then denoise with default
parameters, and remove chimeras:

```sh
vsearch \
    --derep_fulllength reads.fasta \
    --sizeout \
    --output derep.fasta

vsearch \
    --cluster_unoise derep.fasta \
    --sizein \
    --sizeout \
    --centroids denoised.fasta

vsearch \
    --uchime3_denovo denoised.fasta \
    --sizein \
    --nonchimeras amplicons.fasta
```

Denoise with a stricter minimum abundance threshold:

```sh
vsearch \
    --cluster_unoise derep.fasta \
    --sizein \
    --sizeout \
    --minsize 2 \
    --centroids denoised.fasta
```


# SEE ALSO

[`vsearch-cluster_fast(1)`](./vsearch-cluster_fast.1.md),
[`vsearch-cluster_size(1)`](./vsearch-cluster_size.1.md),
[`vsearch-cluster_smallmem(1)`](./vsearch-cluster_smallmem.1.md),
[`vsearch-derep_fulllength(1)`](./vsearch-derep_fulllength.1.md),
[`vsearch-fastx_uniques(1)`](./vsearch-fastx_uniques.1.md),
[`vsearch-uchime3_denovo(1)`](./vsearch-uchime3_denovo.1.md),
[`vsearch-cigar(5)`](../formats/vsearch-cigar.5.md),
[`vsearch-fasta(5)`](../formats/vsearch-fasta.5.md),
[`vsearch-userfields(7)`](../misc/vsearch-userfields.7.md)


#(./fragments/footer.md)
