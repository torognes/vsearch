% vsearch-fastx_getseqs(1) version 2.31.0 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./fragments/date.md)

# NAME

vsearch \-\-fastx_getseqs --- extract sequences from a fasta or fastq
file by label


# SYNOPSIS

| **vsearch** **\-\-fastx_getseqs** _fastxfile_ (**\-\-fastaout** | **\-\-fastqout** | **\-\-notmatched** | **\-\-notmatchedfq**) _filename_ (**\-\-label** _string_ | **\-\-labels** _filename_ | **\-\-label_word** _string_ | **\-\-label_words** _filename_) \[_options_]


# DESCRIPTION

The vsearch command `--fastx_getseqs` extracts sequences from a fasta
or fastq file whose headers match one or more specified labels. The
input format is detected automatically (see
[`vsearch-fasta(5)`](../formats/vsearch-fasta.5.md) and
[`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md)).

This command extends `--fastx_getseq` (see
[`vsearch-fastx_getseq(1)`](./vsearch-fastx_getseq.1.md)) with
additional ways to specify what to match:

- **`--label`** / **`--labels`**: match by full header (or substring
  with `--label_substr_match`). Matching is not case-sensitive.
- **`--label_word`** / **`--label_words`**: match by *word* — a
  sequence of alphanumeric characters delimited by a non-alphanumeric
  character or the start/end of the header. Word matching is
  case-sensitive.
- **`--label_field`**: limit word matching to a named
  semicolon-delimited field in the header (e.g. `abc=123`).

By default, headers are truncated at the first space or tab before
matching — use `--notrunclabels` to disable truncation.

Matching sequences are written to `--fastaout` and/or `--fastqout`.
Non-matching sequences can be written to `--notmatched` and/or
`--notmatchedfq`.

```text
Input (3 sequences):  labels.txt: "seq1"     --labels labels.txt:
                                  "seq3"
>seq1                 extracted (--fastaout):   not matched (--notmatched):
AAAA                  >seq1  AAAA              >seq2  CCCC
>seq2         -->     >seq3  TTTT
CCCC
>seq3
TTTT
```


# OPTIONS

## mandatory options

At least one of `--label`, `--labels`, `--label_word`, or
`--label_words` must be specified. At least one of `--fastaout`,
`--fastqout`, `--notmatched`, or `--notmatchedfq` must also be
specified.


## core options

#(./fragments/option_label.md)

#(./fragments/option_label_field.md)

#(./fragments/option_label_substr_match.md)

#(./fragments/option_label_word.md)

#(./fragments/option_label_words.md)

#(./fragments/option_labels.md)

`--fastaout` *filename*
: Write the matching sequences to *filename*, in fasta format.

`--fastqout` *filename*
: Write the matching sequences to *filename*, in fastq format (see
  [`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md)). Requires fastq
  input.

#(./fragments/option_notmatched.md)

#(./fragments/option_notmatchedfq.md)


## secondary options

#(./fragments/option_bzip2_decompress.md)

#(./fragments/option_fasta_width.md)

#(./fragments/option_fastq_ascii.md)

#(./fragments/option_fastq_qmax.md)

#(./fragments/option_fastq_qmin.md)

#(./fragments/option_gzip_decompress.md)

#(./fragments/option_label_suffix.md)

#(./fragments/option_lengthout.md)

#(./fragments/option_log.md)

#(./fragments/option_no_progress.md)

#(./fragments/option_notrunclabels.md)

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

Extract sequences whose full header matches any label in a list file:

```sh
vsearch \
    --fastx_getseqs input.fasta \
    --labels targets.txt \
    --fastaout extracted.fasta \
    --notmatched rest.fasta
```

Extract sequences by word match, useful when headers contain
annotations after the identifier:

```sh
vsearch \
    --fastx_getseqs input.fasta \
    --label_words targets.txt \
    --fastaout extracted.fasta
```

Extract sequences matching a value in a specific header field:

```sh
vsearch \
    --fastx_getseqs input.fasta \
    --label_field "sample" \
    --label_words samples.txt \
    --fastaout extracted.fasta
```


# SEE ALSO

[`vsearch-fastx_getseq(1)`](./vsearch-fastx_getseq.1.md),
[`vsearch-fastx_getsubseq(1)`](./vsearch-fastx_getsubseq.1.md),
[`vsearch-fasta(5)`](../formats/vsearch-fasta.5.md),
[`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md)


#(./fragments/footer.md)
