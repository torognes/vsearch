% vsearch-chimeras_denovo(1) version 2.30.0 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./fragments/date.md)

# NAME

vsearch \-\-chimeras_denovo --- detect chimeras *de novo* in long exact sequences


# SYNOPSIS

| **vsearch** **\-\-chimeras_denovo** _fastafile_ \[_options_]


# DESCRIPTION

The vsearch command `--chimeras_denovo` detect chimeras *de novo*
(i.e. without external references) in long exact sequences. It uses a
modified uchime algorithm that can automatically adapt to a wide range
of sequence lengths.


# OPTIONS

## mandatory options

At least one of `--alnout`, `--chimeras`, `--nonchimeras`, and
`--tabbedout` must be specified.

`--alnout` *filename*
: Write multi-way alignments to *filename* using a human-readable
  format. Use `--alignwidth` to set the alignment width (60
  nucleotides by default).

#(./fragments/option_chimeras.md)

#(./fragments/option_nonchimeras.md)

`--tabbedout` *filename*
: Write the results to a eighteen-column tab-delimited file with the
  specified *filename*. Columns are:

    1.  score: dummy value, always set to 99.9999
    2.  query header
    3.  parent A header
    4.  parent B header
    5.  parent C header ("*" if there are only two parents)
    6.  QModel: max similarity percentage (always 100.0%)
    7.  QA: similarity percentage with parent A
    8.  QB: similarity percentage with parent B
    9.  QC: similarity percentage with parent C (0.00 if there are only two parents)
    10.  QT: highest similarity percentage with a parent
    11.  left yes: ignored, always set to zero
    12.  left no: ignored, always set to zero
    13.  left abstain: ignored, always set to zero
    14.  right yes: ignored, always set to zero
    15.  right no: ignored, always set to zero
    16.  right abstain: ignored, always set to zero
    17.  dummy value, always set to 0.00
    18.  chimeric status (Y: chimera, N: not a chimera)


## core options

#(./fragments/option_abskew_1.md)

#(./fragments/option_chimeras_diff_pct.md)

#(./fragments/option_chimeras_length_min.md)

#(./fragments/option_chimeras_parents_max.md)

#(./fragments/option_chimeras_parts.md)

#(./fragments/option_sizein.md)

#(./fragments/option_xn.md)


## secondary options

#(./fragments/option_alignwidth.md)

#(./fragments/option_fasta_width.md)

#(./fragments/option_hardmask.md)

#(./fragments/option_label_suffix.md)

#(./fragments/option_log.md)

#(./fragments/option_maxseqlength.md)

#(./fragments/option_minseqlength_32.md)

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

#(./fragments/option_xsize.md)


## pairwise alignment options

These options modify the parameters of the pairwise alignment
model. See
[`vsearch-pairwise_alignment_parameters(7)`](./misc/vsearch-pairwise_alignment_parameters.7.md)
for a complete description. Modify with caution.

#(./fragments/option_gapext.md)

#(./fragments/option_gapopen.md)

#(./fragments/option_match.md)

#(./fragments/option_mismatch.md)


## ignored options

#(./fragments/option_threads_not_multithreaded.md)


## unsupported options

The following options are not yet supported by the `--chimeras_denovo`
command:

#(./fragments/option_bzip2_decompress.md)

#(./fragments/option_gzip_decompress.md)

#(./fragments/option_lengthout.md)

#(./fragments/option_xlength.md)


# EXAMPLES


# SEE ALSO

[`vsearch-uchime_denovo(1)`](./commands/vsearch-uchime_denovo.1.md),
[`vsearch-uchime2_denovo(1)`](../formats/vsearch-uchime2_denovo.1.md),
[`vsearch-uchime3_denovo(1)`](../formats/vsearch-uchime3_denovo.1.md),
[`vsearch-uchime_ref(1)`](../formats/vsearch-uchime_ref.1.md)


#(./fragments/footer.md)
