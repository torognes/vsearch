% vsearch-cut(1) version 2.30.0 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./fragments/date.md)

# NAME

vsearch \-\-cut --- use a restriction pattern to cut fasta sequences


# SYNOPSIS

| **vsearch** **\-\-cut** *fastafile* \-\-cut_pattern *pattern* (\-\-fastaout | \-\-fastaout_rev | \-\-fastaout_discarded | \-\-fastaout_discarded_rev) *outputfile* \[*options*]


# DESCRIPTION

The vsearch command `--cut` uses a restriction pattern to cut input
fasta sequences. Input sequences are cut into fragments at **each**
restriction site matching the pattern given with the option
`--cut_pattern`. Restriction patterns are only searched on the forward
(or normal) strand, not on the reverse strand.

Fragments on the forward strand are written to the file specified with
the `--fastaout` file, and the reverse-complement of fragments are
written to the file specified with the `--fastaout_rev`
option. Fragments receive the name of their parent sequence. Input
sequences with not match are written to the file specified with the
option `--fastaout_discarded`, and their reverse-complement are also
written to the file specified with the `--fastaout_discarded_rev`
option.

A typical restriction pattern is "`G^AATT_C`", representing the EcoRI
restriction site. The nucleotide symbols represent the sequence to be
matched. Lowercase or uppercase nucleotides, as well as ambiguous
nucleotides (IUPAC) are accepted. The special character '^'
(circumflex) indicates the cutting position on the forward strand,
while '`_`' (underscore) indicates the cutting position on the reverse
strand. Forward and reverse cutting positions can be the same (for
example "`GG^_AA`"), but exactly one cutting position on each strand
must be indicated (one '^' and one '`_`').

As noted above, restriction patterns are only searched on the forward
(or normal) strand, not on the reverse strand. For palindromic
patterns such as EcoRI, this is not an issue. A palindromic pattern is
identical to its reverse-complement (same sequence and same cutting
positions), so there are no copies on the reverse strand that do not
have a counterpart on the forward strand.

```
forward 5' A...G^AATT_C...G 3'
reverse 3' T...C_TTAA^G...C 5'

fragments: A...G, AATTC...G, but also C...G and AATTC...T
```

All copies are detected.

For asymetrical or non-palindromic patterns, if the pattern appears on
the reverse strand, it will not be detected. For example, the pattern
"`GG^_C`" is detected when present on the forward strand:

```
forward 5' A...GG^C...G 3'
reverse 3' T...CC G...C 5'

fragments: A...GG, C...G
```

but is not detected when present on the reverse strand:

```
forward 5' A...G CC...G 3'
reverse 3' T...C_GG...C 5'  <- present but not detected
```

To detect asymetrical or non-palindromic patterns on the reverse
strand, it is necessary to run the `--cut` command on the
reverse-complemented input sequences (see
[`vsearch-fastx_revcomp(1)`](./vsearch-fastx_revcomp.1.md)).

Finally, pattern occurrences can overlap. For example the pattern
"`G^_G`" will be detected *three* times in the sequence "`GGGG`".


# OPTIONS

## mandatory options

`--cut_pattern` *pattern*
: Specify the restriction site pattern (case insensitive IUPAC
  characters) and cutting positions ('^' and '_'). For example,
  "`G^AATT_C`". See the PATTERN EXAMPLES section for more details.

At least one of `--fastaout`, `--fastaout_rev`,
`--fastaout_discarded`, or `--fastaout_discarded_rev` must also be
specified.


## core options

`--fastaout` *filename*
: Write the forward strand fragments to *filename*, in fasta format.

`--fastaout_rev` *filename*
: Write the reverse strand fragments to *filename*, in fasta format.

`--fastaout_discarded` *filename*
: Write the non-matching sequences to *filename*, in fasta format.

`--fastaout_discarded_rev` *filename*
: Write the reverse-complemented non-matching sequences to *filename*,
  in fasta format.


## secondary options

#(./fragments/option_bzip2_decompress.md)

#(./fragments/option_fasta_width.md)

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

Cut the sequences in *query.fasta*, using the restriction pattern
"`G^AATT_C`" (`--cut_pattern`). Write the fragments found on the
reverse strand to *query_fragments.fasta*, in fasta format
(`--fastaout_rev`):

```sh
vsearch \
    --cut query.fasta \
    --cut_pattern "G^AATT_C" \
    --fastaout_rev query_fragments.fasta
```


# PATTERN EXAMPLES

(the symbol '|' is used to represent a cut on either strand)

`EcoRI`
: use the palindromic pattern "`G^AATT_C`" to represent the following
  restriction:

```
5' G|AATT-C 3'
3' C-TTAA|G 5'
```

`EcoRII`
: use the palindromic pattern "`^CCWGG_`" to represent the following
  restriction:

```
5' |CCWGG- 3'
3' -GGWCC| 5'
```

`EcoRV`
: use the palindromic pattern "`GAT^_ATC`" to represent the following
  restriction:

```
5' GAT|ATC 3'
3' CTA|TAG 5'
```

`Fok1`
: use the non-palindromic pattern "`N_NNNNNNNNNNNNGGATGNNNNNNNN^N`" to
represent the following restriction, with a cut 9 nucleotides
downstream of the motif on the forward strand, and 13 nucleotides
downstream of the motif on the reverse strand:

```
5' -NNNNNNNNNNNNGGATGNNNNNNNN| 3'
3' |NNNNNNNNNNNNCCTACNNNNNNNN- 5'
```

`HindIII`
: use the palindromic pattern "`A^AGCT_T`" to represent the following
  restriction:

```
5' A|AGCT-T 3'
3' T-TCGA|A 5'
```

`NspI`
: use the asymetric pattern "`RCATG^_Y`" to represent the following
  restriction:

```
5' RCATG|Y 3'
3' YGTAC-R 5'
```


# SEE ALSO

[`vsearch-fastx_revcomp(1)`](./vsearch-fastx_revcomp.1.md),


#(./fragments/footer.md)
