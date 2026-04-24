% vsearch-uchime3_denovo(1) version 2.31.0 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./fragments/date.md)

# NAME

vsearch \-\-uchime3_denovo --- detect chimeras *de novo* using the UCHIME2 algorithm with a stricter abundance skew


# SYNOPSIS

| **vsearch** **\-\-uchime3_denovo** _fastafile_ (**\-\-chimeras** | **\-\-nonchimeras** | **\-\-uchimealns** | **\-\-uchimeout**) _filename_ \[_options_]


# DESCRIPTION

The vsearch command `--uchime3_denovo` detects chimeric sequences
present in the fasta-formatted *fastafile*, without the use of
external references (*de novo*). It uses the UCHIME2 algorithm, which
is designed for denoised amplicons (see `--cluster_unoise`). The only
difference from `--uchime2_denovo` is that the default minimum
abundance skew (`--abskew`) is set to 16.0 rather than 2.0, reflecting
the higher abundance contrast expected between genuine sequences and
chimeras in well-denoised datasets. Sequences are compared on their
*plus* strand only.

Chimera detection is based on a scoring function controlled by two
options: `--dn` and `--xn`. Note that `--mindiffs`, `--mindiv`, and
`--minh` are ignored by this command.

Input sequences must carry abundance annotations in their headers
(e.g. `;size=integer;`). `--sizein` is always implied; it does not
need to be specified. Sequences are automatically sorted by decreasing
abundance before chimera detection. The assumption is that chimeras
appear later in the PCR amplification process and are therefore less
abundant than their parents (see `--abskew`).

At least one output option must be specified.

See also `--uchime_denovo` for the original UCHIME algorithm and
`--uchime_ref` for reference-based chimera detection.


# OPTIONS

## mandatory options

`--uchime3_denovo` *fastafile*
: Detect chimeras *de novo* in the fasta-formatted *fastafile* using
  the UCHIME2 algorithm with a stricter default abundance skew. This
  option is mandatory.

At least one of the following output options must be specified:

#(./fragments/option_chimeras.md)

#(./fragments/option_nonchimeras.md)

#(./fragments/option_uchimealns.md)

#(./fragments/option_uchimeout.md)

The `--borderline` option can also produce output, but only when
combined with at least one of the output options listed above:

#(./fragments/option_borderline.md)


## core options

#(./fragments/option_abskew_16.md)

#(./fragments/option_dn.md)

#(./fragments/option_sizein.md)
: Always implied.

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

#(./fragments/option_mindiffs.md)
: Ignored by `--uchime3_denovo`.

#(./fragments/option_mindiv.md)
: Ignored by `--uchime3_denovo`.

#(./fragments/option_minh.md)
: Ignored by `--uchime3_denovo`.

#(./fragments/option_threads_not_multithreaded.md)


# EXAMPLES

Detect chimeras *de novo* using the UCHIME2 algorithm with a strict
abundance skew, and write non-chimeric sequences to a file:

```sh
vsearch \
    --uchime3_denovo denoised.fasta \
    --sizein \
    --nonchimeras clean.fasta
```

A typical denoising workflow: dereplicate, denoise, then remove
chimeras (recommended approach after `--cluster_unoise`):

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


# SEE ALSO

[`vsearch-uchime_denovo(1)`](./vsearch-uchime_denovo.1.md),
[`vsearch-uchime2_denovo(1)`](./vsearch-uchime2_denovo.1.md),
[`vsearch-uchime_ref(1)`](./vsearch-uchime_ref.1.md),
[`vsearch-chimeras_denovo(1)`](./vsearch-chimeras_denovo.1.md),
[`vsearch-cluster_unoise(1)`](./vsearch-cluster_unoise.1.md),
[`vsearch-derep_fulllength(1)`](./vsearch-derep_fulllength.1.md),
[`vsearch-sortbysize(1)`](./vsearch-sortbysize.1.md),
[`vsearch-fasta(5)`](../formats/vsearch-fasta.5.md)


#(./fragments/footer.md)
