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
`ascii(7)`). The control characters 1–8, 11–12, 14–31 and 127 cause a
fatal error; the tabulation (9) is accepted; a NUL byte (0) silently
truncates the label; carriage returns (13) and newlines (10) terminate
the line. Non-ASCII characters (values 128–255) trigger a warning.
This check covers the part of the header the command retains: the
label by default, or the whole header line for the filtering commands
(`--fastq_filter`, `--fastx_filter`) and when `--notrunclabels` is in
effect.

## Sequence

#(./fragments/format_sequence_fastq.md)

## Quality header

The *quality header* is the line starting with '+'. The text after
'+' must be empty or repeat the sequence header line exactly; anything
else causes a fatal error.

## Quality string

The *quality string* is a string of ASCII characters, one per base,
encoding the Phred quality score for the corresponding position. The
quality string starts after the quality header line and ends before the
next sequence header or the end of the file.

Quality scores are encoded as `Q = ASCII_value - offset`, where the
offset is set with `--fastq_ascii` (default 33):

- **phred+33** (offset 33): Sanger and Illumina 1.8+ format. Valid
  quality characters range from '!' (Q=0) to '~' (Q=93). vsearch
  accepts the whole range by default (`--fastq_qmin` 0 and
  `--fastq_qmax` 93). Before version 3.0 the upper bound defaulted to
  41, which rejected PacBio HiFi and nanopore files.

- **phred+64** (offset 64): Illumina 1.3+ and Illumina 1.5+ formats.
  Valid quality characters range from '@' (Q=0) to '~' (Q=62).

Because the accepted range now covers every representable score, an
out-of-range error no longer doubles as a wrong-`--fastq_ascii` detector:
a phred+64 file read at offset 33 decodes to legal Sanger scores 31 too
high, and would once have been stopped by the old bound of 41. vsearch
therefore warns when the quality symbols it read contradict the requested
offset, using the same heuristic `--fastq_chars` prints its guess from.
The warning is advisory: the two encodings overlap, since a file whose
scores all exceed Q30 is indistinguishable from a phred+64 file by its
symbols alone. It is raised only for files of at least 100 records, below
which the observed range is not evidence, and `--fastq_chars` does not
raise it because it reports its own guess. Run `vsearch --fastq_chars
FILE` whenever it appears.

The older Solexa/Illumina 1.0 format shares the offset 64 but not the
score definition, and **vsearch does not support it**. A Solexa score is
`Q = -10 log10(p / (1 - p))`, not the Phred `Q = -10 log10(p)`, and in
the Illumina 1.0 pipeline it ranges from -5 to 40 (';' to 'h'). vsearch
always applies the Phred formula, so a Solexa quality string is read as
if it were a Phred one: the error probability comes out 1.3 times too
large at Q=5, 2 times at Q=0, and 4.2 times at Q=-5, where the Phred
formula returns 3.16 --- not a probability at all. The two definitions
converge as the score rises (3% apart at Q=15, 1% at Q=20, 0.1% at
Q=30), so only the low scores are affected, but those are the ones a
quality filter acts on. Convert such files to phred+33 with a dedicated
tool before feeding them to vsearch.

The characters themselves are not rejected. Scores below zero are
representable down to -31 by lowering `--fastq_qmin` (the offset plus
the minimum must stay at 33 or more), and `--fastq_chars` diagnoses a
Solexa file correctly. What the commands then do with a negative score
differs: `--fastq_convert` clamps it to `--fastq_qminout` (0 by
default); `--fastq_eestats` and `--fastq_eestats2` raise it to 0;
`--fastq_mergepairs` and `--fastx_uniques` treat every score below 2 as
a blind guess (p = 0.75), and `--fastx_uniques` writes it out as Q=1;
and `--fastx_filter` discards the whole sequence, since
`--fastq_minqual` defaults to 0 and may not be negative.

Which encoding an *output* file carries depends on the command.
`--fastq_convert`, `--fastx_uniques`, `--fasta2fastq` and `--sff_convert`
re-encode, using `--fastq_asciiout` (default 33) independently of the
input offset, and clamping to `--fastq_qminout` and `--fastq_qmaxout`.
`--fastq_mergepairs` writes the qualities it computes with the *input*
offset, `--fastq_ascii`, and does not accept `--fastq_asciiout`. Every
other command copies the quality symbols verbatim, so the output carries
whatever encoding the input had and `--fastq_asciiout` is not accepted
either. To change the encoding of a file, use `--fastq_convert`.

The parser accepts every printable ASCII character (values 33–126) in
quality strings --- '.' and '-' are ordinary quality symbols (Q=13 and
Q=12 at offset 33). Carriage returns before the newline are accepted;
any other character (tabulations, spaces, control characters, values
128–255) causes a fatal error; nothing is stripped or warned about.
The `--fastq_qmin`/`--fastq_qmax` score bounds are enforced by the
commands that interpret quality values, not by the parser itself.

## Header annotations

vsearch reads and writes several annotations embedded in sequence headers,
following the pattern `[@;]key=value[;]`. Annotations can appear in any
order, but must be joined to the label by ';' without any whitespace:
with the default header truncation, an annotation separated from the
label by a space is part of the discarded remainder and is silently
invisible (unless `--notrunclabels` is used).

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
