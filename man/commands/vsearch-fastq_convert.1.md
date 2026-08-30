% vsearch-fastq_convert(1) version 2.31.0 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./fragments/date.md)

# NAME

vsearch \-\-fastq_convert --- convert between fastq format variants


# SYNOPSIS

| **vsearch** **\-\-fastq_convert** _fastqfile_ **\-\-fastqout** _filename_ \[_options_]


# DESCRIPTION

The vsearch command `--fastq_convert` converts a fastq file from one
quality-score encoding to another. The input encoding is specified with
`--fastq_ascii` (default: 33) and the output encoding with
`--fastq_asciiout` (default: 33). Both accept the values 33 (phred+33,
Sanger/Illumina 1.8+) and 64 (phred+64, Illumina 1.3+/1.5+). The older
Solexa/Illumina 1.0 format shares the offset 64 but defines its scores
differently (see
[`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md)); `--fastq_solexa`
converts such a file to the Phred scale, and is the only place in
vsearch where the Solexa score definition is understood.

Quality scores are remapped during conversion. Output scores are
always clamped to the range set by `--fastq_qminout` and
`--fastq_qmaxout`, whose defaults are 0 and 93: scores are written
unchanged unless `--fastq_qmaxout` is lowered (before version 3.0 the
default was 41, and any higher score was silently reduced to it), and
negative scores (e.g. Solexa scores read
with `--fastq_ascii 64 --fastq_qmin -5`) are raised to 0, unless the
bounds are changed. The input score range is validated against
`--fastq_qmin` and `--fastq_qmax`.

Use `--fastq_chars` (see
[`vsearch-fastq_chars(1)`](./vsearch-fastq_chars.1.md)) to detect the
encoding of an unknown fastq file before converting.

To illustrate a conversion from phred+64 to phred+33:

```text
Input (phred+64):    Output (--fastq_ascii 64 --fastq_asciiout 33):

@s1                  @s1
ACGT                 ACGT
+                    +
hijk                 IJKL    (same Phred scores, different ASCII offset)
```

## Converting Solexa scores

A Solexa score is *Q = -10 log10(p / (1 - p))*, while every score
vsearch understands is a Phred score, *Q = -10 log10(p)*. The two share
the ASCII offset 64 and nothing else, so rebasing the offset alone
leaves the scores wrong. `--fastq_solexa` applies

```text
Q_phred = 10 log10(10^(Q_solexa / 10) + 1)
```

and rounds the result to the nearest integer. It requires
`--fastq_ascii 64`, and implies `--fastq_qmin -5` (the floor of the
Solexa scale) unless `--fastq_qmin` is given explicitly.

The mapping is the identity from Solexa 10 upward, so only the fifteen
lowest scores change --- which are exactly the ones a quality filter
acts on:

```text
Solexa:  -5  -4  -3  -2  -1   0   1   2   3   4   5   6   7   8   9  10
Phred:    1   1   2   2   3   3   4   4   5   5   6   7   8   9  10  10
```

The conversion is one-way and lossy: six Solexa pairs collapse onto a
single Phred score each ({-5, -4}, {-3, -2}, {-1, 0}, {1, 2}, {3, 4}
and {9, 10}), so the original scores cannot be recovered afterwards.
There is no option to write Solexa scores. This is inherent to the two
scales rather than to the implementation.


# OPTIONS

## mandatory options

#(./fragments/option_fastqout_convert.md)


## core options

#(./fragments/option_fastq_ascii.md)

#(./fragments/option_fastq_asciiout.md)

#(./fragments/option_fastq_qmaxout.md)

#(./fragments/option_fastq_qminout.md)

#(./fragments/option_fastq_solexa.md)


## secondary options

#(./fragments/option_bzip2_decompress.md)

#(./fragments/option_fastq_qmax.md)

#(./fragments/option_fastq_qmin.md)

#(./fragments/option_gzip_decompress.md)

#(./fragments/option_label_suffix.md)

#(./fragments/option_lengthout.md)

#(./fragments/option_log.md)

#(./fragments/option_no_progress.md)

#(./fragments/option_quiet.md)

#(./fragments/option_relabel.md)

#(./fragments/option_relabel_keep.md)

#(./fragments/option_relabel_md5.md)

#(./fragments/option_relabel_self.md)

#(./fragments/option_relabel_sha1.md)

#(./fragments/option_sample.md)

#(./fragments/option_sizein.md)

#(./fragments/option_sizeout.md)

#(./fragments/option_xee.md)

#(./fragments/option_xlength.md)

#(./fragments/option_xsize.md)


## ignored options

#(./fragments/option_threads_not_multithreaded.md)


# EXAMPLES

Convert a phred+64 file to the standard phred+33 encoding:

```sh
vsearch \
    --fastq_convert input.fastq \
    --fastq_ascii 64 \
    --fastq_asciiout 33 \
    --fastqout converted.fastq
```

Convert and clamp output quality scores to the range 0--40 (an older
Illumina convention):

```sh
vsearch \
    --fastq_convert input.fastq \
    --fastq_ascii 64 \
    --fastq_asciiout 33 \
    --fastq_qmaxout 40 \
    --fastqout converted.fastq
```

Convert a Solexa/Illumina 1.0 file to phred+33, converting the score
definition as well as the offset:

```sh
vsearch \
    --fastq_convert input.fastq \
    --fastq_ascii 64 \
    --fastq_solexa \
    --fastqout converted.fastq
```


# SEE ALSO

[`vsearch-fastq_chars(1)`](./vsearch-fastq_chars.1.md),
[`vsearch-fastq_stats(1)`](./vsearch-fastq_stats.1.md),
[`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md)


#(./fragments/footer.md)
