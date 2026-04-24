% vsearch-fastq(5) version 2.31.0 | vsearch file formats
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(../commands/fragments/date.md)

# NAME

fastq --- a text-based format for representing nucleotide sequences
and their corresponding quality scores


# DESCRIPTION

In fastq files, each entry consists of four parts: a *sequence header*
starting with '@', a *sequence*, a *quality header* starting with '+',
and a *quality string*:

```text
@label [description]
ACGTURYSWKMDBHVN...
+
IIIIIIIIIIIIIIII...
```

## Sequence header

The *sequence header* is the line starting with '@'. The *label* is the
string between '@' and the first space, tab, or newline, unless
`--notrunclabels` is in effect, in which case the entire line after '@'
is used as the label.

The header should contain printable ASCII characters (values 33–126, see
`ascii(7)`). Unprintable characters (values 0–32 and 127) cause a fatal
error. Non-ASCII characters (values 128–255) trigger a warning.

## Sequence

#(./fragments/format_sequence.md)

## Quality header

The *quality header* is the line starting with '+'. Any text after '+'
on this line is ignored by vsearch.

## Quality string

The *quality string* is a string of ASCII characters, one per base,
encoding the Phred quality score for the corresponding position. The
quality string starts after the quality header line and ends before the
next sequence header or the end of the file.

Quality scores are encoded as `Q = ASCII_value - offset`, where the
offset is set with `--fastq_ascii` (default 33):

- **phred+33** (offset 33): Sanger and Illumina 1.8+ format. Valid
  quality characters range from '!' (Q=0) to '~' (Q=93). vsearch
  accepts scores from 0 to 41 by default (`--fastq_qmin` and
  `--fastq_qmax`).

- **phred+64** (offset 64): Solexa, Illumina 1.3+, and Illumina 1.5+
  formats. Valid quality characters range from '@' (Q=0) to '~' (Q=62).

vsearch silently ignores ASCII characters 9–13 (tab and whitespace) in
quality strings, and exits with an error if ASCII characters 0–8, 14–31,
'.' or '-' are present. All other unexpected characters are stripped with
a warning.

## Header annotations

vsearch reads and writes several annotations embedded in sequence headers,
following the pattern `[@;]key=value[;]`. Annotations can appear in any
order after the label:

`[@;]size=integer[;]`
: Abundance (number of occurrences of the sequence in the study). Read by
  `--sizein`, written by `--sizeout`, removed by `--xsize`.

`[@;]ee=float[;]`
: Expected error count. Written by `--eeout` or `--fastq_eeout`, removed
  by `--xee`. See
  [`vsearch-expected_error(7)`](../misc/vsearch-expected_error.7.md).

`[@;]length=integer[;]`
: Sequence length. Written by `--lengthout`, removed by `--xlength`.

`[@;]sample=string[;]`
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

[`vsearch-fasta(5)`](./vsearch-fasta.5.md),
[`vsearch-udb(5)`](./vsearch-udb.5.md),
[`vsearch-fastq_chars(1)`](../commands/vsearch-fastq_chars.1.md),
[`vsearch-fastq_convert(1)`](../commands/vsearch-fastq_convert.1.md),
`ascii(7)`


#(../commands/fragments/footer.md)
