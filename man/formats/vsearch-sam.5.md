% vsearch-sam(5) version 2.31.0 | vsearch file formats
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(../commands/fragments/date.md)

# NAME

sam --- a text-based format for representing sequence alignments


# DESCRIPTION

SAM (Sequence Alignment/Map) is a tab-delimited text format for
storing the alignment of one or more *query* sequences against one or
more *reference* sequences. It was designed by the 1000 Genomes
Project and is maintained by the hts-specs working group:

<https://github.com/samtools/hts-specs>

The format is defined in the *SAMv1* specification, and optional
alignment tags are listed in the *SAMtags* specification. vsearch
produces SAM version 1.0 and writes each record on a single line with
fields separated by a single horizontal tab character (ASCII 9). Lines
are terminated by a single line feed character (ASCII 10).

A SAM file produced by vsearch consists of an optional *header*
followed by zero or more *alignment records*. The `--samout` *filename*
option selects SAM as the output format for the `--allpairs_global`,
`--cluster_fast`, `--cluster_size`, `--cluster_smallmem`,
`--cluster_unoise`, `--search_exact`, and `--usearch_global` commands.
The `--samheader` option requests that the optional header be written;
by default no header is produced.

vsearch does not read SAM files and does not produce the binary BAM or
CRAM formats. Conversion between SAM, BAM, and CRAM can be performed
with external tools such as `samtools(1)`.


## Header

When `--samheader` is given, vsearch writes a header consisting of one
`@HD` line, one `@SQ` line per reference sequence, and one `@PG` line,
in that order. All header lines start with the character `@`
immediately followed by a two-letter record type code, and contain
tab-separated `TAG:value` fields:

```text
@HD	VN:1.0	SO:unsorted	GO:query
@SQ	SN:name	LN:length	M5:md5	UR:file:dbname
@PG	ID:vsearch	VN:version	CL:command-line
```

`@HD`
: File-level metadata. `VN` is the SAM format version (`1.0`). `SO` is
  the sort order of alignment records (`unsorted`). `GO` is the
  grouping of alignment records (`query`: records for the same query
  are consecutive).

`@SQ`
: Reference sequence dictionary. One `@SQ` line is emitted per
  sequence in the reference database, in database order. `SN` is the
  reference name, taken from the fasta header up to the first white
  space (or the entire header line with `--notrunclabels`). `LN` is
  the reference length in residues. `M5` is the MD5 digest of the
  reference sequence, as 32 lowercase hexadecimal digits. `UR` is the
  location of the reference, written as `file:` followed by the
  database filename as given on the command line.

`@PG`
: Program record. `ID` is the program name (`vsearch`). `VN` is the
  vsearch version. `CL` is the full command line used to invoke
  vsearch.

vsearch does not emit `@RG` (read group) or `@CO` (comment) header
lines.


## Alignment records

Each alignment record occupies a single line with eleven mandatory
tab-separated fields followed by zero or more optional tags:

```text
QNAME	FLAG	RNAME	POS	MAPQ	CIGAR	RNEXT	PNEXT	TLEN	SEQ	QUAL	[TAG:TYPE:VALUE ...]
```

1. `QNAME` — query template name. The query label as it appears in
   the input fasta or fastq file.
2. `FLAG` — bitwise flag, written as a decimal integer. vsearch only
   sets the following bits:
   - `0x004` (`4`): the query is unmapped (no alignment was produced).
     Set when the query has no hit and `--output_no_hits` is in
     effect.
   - `0x010` (`16`): the query aligned on the reverse strand of the
     reference. Set for reverse-strand hits produced with `--strand
     both`.
   - `0x100` (`256`): the record is a secondary alignment. Set for the
     second and following hits reported for a given query when
     `--maxaccepts` and `--top_hits_only` allow multiple hits.
   All other bits are always `0`; in particular, vsearch never sets
   the paired-end bits (`0x001`, `0x002`, `0x040`, `0x080`), the mate
   bits (`0x008`, `0x020`), the quality-control bit (`0x200`), the
   duplicate bit (`0x400`), or the supplementary bit (`0x800`).
3. `RNAME` — reference sequence name. The target label as it appears
   in the reference database, matching the `SN` value of the
   corresponding `@SQ` header line. `*` for unmapped queries.
4. `POS` — 1-based leftmost mapping position of the query on the
   reference. Always `1` for mapped queries (vsearch performs global
   alignment and reports the alignment relative to the start of the
   reference). `0` for unmapped queries.
5. `MAPQ` — mapping quality. Always `255`, the value defined by the
   SAM specification to mean *mapping quality not available*.
6. `CIGAR` — alignment operations in CIGAR format, with the operation
   letters `M` (match or mismatch), `I` (insertion in the query
   relative to the reference), and `D` (deletion in the query relative
   to the reference). vsearch always emits explicit run-lengths and
   does not use the `=` and `X` operations. `*` for unmapped queries.
   See [`vsearch-cigar(5)`](./vsearch-cigar.5.md).
7. `RNEXT` — reference sequence name of the mate. Always `*` (vsearch
   does not output paired-end information).
8. `PNEXT` — position of the mate. Always `0`.
9. `TLEN` — observed template length. Always `0`.
10. `SEQ` — query segment. For forward-strand hits, the original query
    sequence. For reverse-strand hits (flag bit `0x010` set), the
    reverse complement of the original query sequence. For unmapped
    queries, the original query sequence. vsearch never uses `*` in
    this field.
11. `QUAL` — ASCII-encoded Phred quality scores (base 33) for each
    residue in `SEQ`. Always `*` (not available), even when the input
    is a fastq file.

For each mapped query, vsearch emits one alignment record per
reported hit. The number and ordering of hits are controlled by
`--maxaccepts`, `--maxrejects`, `--top_hits_only`, and other
search-related options. The first record of a group has flag bit
`0x100` cleared and represents the primary alignment; subsequent
records for the same query have flag bit `0x100` set and represent
secondary alignments.


## Optional tags

After the eleven mandatory fields, vsearch appends the following
optional tags, in order, for each mapped query. Each tag is written
as `TAG:TYPE:VALUE` and separated from the previous field by a tab.
Unmapped queries carry no optional tags.

`AS:i:`*integer*
: Alignment score. vsearch stores the percent identity of the
  alignment, rounded to the nearest integer. This differs from the
  SAMtags recommendation, which describes `AS` as an aligner-specific
  alignment score; vsearch follows the convention established by
  usearch.

`XN:i:0`
: Next best alignment score. Always `0` (not computed by vsearch).

`XM:i:`*integer*
: Number of mismatches in the alignment, excluding terminal gaps.

`XO:i:`*integer*
: Number of gap opens in the alignment, excluding terminal gaps.

`XG:i:`*integer*
: Total number of gap extensions in the alignment, excluding terminal
  gaps. Equivalent to the total length of internal gaps.

`NM:i:`*integer*
: Edit distance to the reference, defined as the sum of `XM` and
  `XG` (mismatches plus total internal gap length).

`MD:Z:`*string*
: Mismatching positions and deleted reference bases, encoded as
  described in the SAMtags specification. The `MD` string allows the
  reference sequence within the aligned region to be reconstructed
  from `SEQ` and `CIGAR`. It consists of a sequence of decimal numbers
  (counts of consecutive matching bases), uppercase letters (mismatched
  reference bases), and `^`-prefixed uppercase letter runs (reference
  bases deleted from the query). Insertions in the query are not
  represented in `MD`; they are recorded only in `CIGAR`.

`YT:Z:UU`
: Alignment type in bowtie2 notation. Always `UU` (unpaired read
  aligned to a single unpaired reference). vsearch does not align
  paired reads.

The `XN`, `XM`, `XO`, `XG`, `MD`, and `YT` tags are not part of the
SAMv1 standard but are widely used and originate from bowtie2 and
usearch.


## Unmapped queries

By default, a query that produces no hit does not generate a SAM
record. When `--output_no_hits` is given, vsearch emits a single
record for each such query with `FLAG` set to `4` (unmapped), `RNAME`
set to `*`, `POS` set to `0`, `CIGAR` set to `*`, and no optional
tags. The `SEQ` field contains the original query sequence.


## Character set and line endings

vsearch writes ASCII characters only. Lines are terminated by a single
line feed character (ASCII 10); no carriage return is emitted. Fields
and tags are separated by a single tab character (ASCII 9). vsearch
does not pad fields and does not emit trailing tabs.


# EXAMPLES

A minimal SAM file with header, one forward-strand hit, one
reverse-strand hit for a second query, and one unmapped query
(`--output_no_hits` in effect):

```text
@HD	VN:1.0	SO:unsorted	GO:query
@SQ	SN:ref1	LN:10	M5:0123456789abcdef0123456789abcdef	UR:file:refs.fasta
@SQ	SN:ref2	LN:12	M5:fedcba9876543210fedcba9876543210	UR:file:refs.fasta
@PG	ID:vsearch	VN:2.30.6	CL:vsearch --usearch_global queries.fasta --db refs.fasta --id 0.9 --samout out.sam --samheader
q1	0	ref1	1	255	10M	*	0	0	ACGTACGTAC	*	AS:i:100	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:10	YT:Z:UU
q2	16	ref2	1	255	5M1I6M	*	0	0	ACGTAGCGTACG	*	AS:i:92	XN:i:0	XM:i:0	XO:i:1	XG:i:1	NM:i:1	MD:Z:11	YT:Z:UU
q3	4	*	0	255	*	*	0	0	TTTTTTTT	*
```

The second record has the reverse-strand bit `0x010` (decimal `16`)
set, meaning that `SEQ` is the reverse complement of the original
query and the alignment lies on the reverse strand of `ref2`. The
third record is unmapped.


# SEE ALSO

[`vsearch-allpairs_global(1)`](../commands/vsearch-allpairs_global.1.md),
[`vsearch-cluster_fast(1)`](../commands/vsearch-cluster_fast.1.md),
[`vsearch-cluster_size(1)`](../commands/vsearch-cluster_size.1.md),
[`vsearch-cluster_smallmem(1)`](../commands/vsearch-cluster_smallmem.1.md),
[`vsearch-cluster_unoise(1)`](../commands/vsearch-cluster_unoise.1.md),
[`vsearch-search_exact(1)`](../commands/vsearch-search_exact.1.md),
[`vsearch-usearch_global(1)`](../commands/vsearch-usearch_global.1.md),
[`vsearch-cigar(5)`](./vsearch-cigar.5.md),
[`vsearch-fasta(5)`](./vsearch-fasta.5.md),
[`vsearch-fastq(5)`](./vsearch-fastq.5.md)

The SAM and SAMtags specifications are maintained at
<https://github.com/samtools/hts-specs>.


#(../commands/fragments/footer.md)
