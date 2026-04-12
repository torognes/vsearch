% vsearch-fasta(5) version 2.30.6 | vsearch file formats
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(../commands/fragments/date.md)

# NAME

fasta --- a text-based format for representing nucleotide sequences


# DESCRIPTION

In fasta files, each entry consists of a *header* line and one or more
*sequence* lines:

```text
>label [description]
ACGTURYSWKMDBHVN...
```

## Header

The *header* is the line starting with '>'. The *label* is the string
between '>' and the first space, tab, or newline, unless `--notrunclabels`
is in effect, in which case the entire line after '>' is used as the label.

The header should contain printable ASCII characters (values 33–126, see
`ascii(7)`). Unprintable characters (values 0–32 and 127) cause a fatal
error. Non-ASCII characters (values 128–255) trigger a warning.

## Sequence

#(./fragments/format_sequence.md)

## Header annotations

vsearch reads and writes several annotations embedded in sequence headers,
following the pattern `[>;]key=value[;]`. Annotations can appear in any
order after the label:

`[>;]size=integer[;]`
: Abundance (number of occurrences of the sequence in the study). Read by
  `--sizein`, written by `--sizeout`, removed by `--xsize`.

`[>;]ee=float[;]`
: Expected error count. Written by `--eeout` or `--fastq_eeout`, removed
  by `--xee`. See
  [`vsearch-expected_error(7)`](../misc/vsearch-expected_error.7.md).

`[>;]length=integer[;]`
: Sequence length. Written by `--lengthout`, removed by `--xlength`.

`[>;]sample=string[;]`
: Sample identifier. Written by `--sample`.

## Case and masking

vsearch operations are case-insensitive, except when soft masking is
active. Lower-case letters indicate masked residues; upper-case letters
indicate unmasked residues. DUST masking (automatic low-complexity masking
applied during chimera detection, clustering, masking, pairwise alignment,
and searching) converts low-complexity regions to lower case. Soft masking
is enabled with `--qmask soft` or `--dbmask soft`. Masked residues are
excluded from word-based comparisons but otherwise treated as normal
symbols.

## T and U equivalence

When comparing sequences (chimera detection, dereplication, searching, and
clustering), T and U are treated as identical regardless of case.


# SEE ALSO

[`vsearch-fastq(5)`](./vsearch-fastq.5.md),
[`vsearch-udb(5)`](./vsearch-udb.5.md)


#(../commands/fragments/footer.md)
