% vsearch-usearch(7) version 2.32.0 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(../commands/fragments/date.md)

# NAME

usearch compatibility --- the commands, options and behaviours vsearch
shares with usearch, and those it does not


# DESCRIPTION

If you are a usearch user, our objective is to make you feel at home:
vsearch was designed to behave like usearch, to some extent. Like any
complex software, usearch is not free from quirks and inconsistencies.
We decided not to reproduce some of them, and, for complete
transparency, this page documents the differences that matter when a
usearch command line is ported to vsearch: commands vsearch does not
implement, options it spells differently, options it accepts and
ignores, options it rejects, and behaviours it changed on purpose.

Two differences are broad enough to state first. vsearch works with
nucleotide sequences only; amino acid sequences are not supported (see
[`vsearch-nucleotides(7)`](./vsearch-nucleotides.7.md)). And vsearch
always aligns globally, with full dynamic programming, so the options
by which usearch tunes or bounds its alignment stage have nothing to
act on (see
[`vsearch-pairwise_alignment_parameters(7)`](./vsearch-pairwise_alignment_parameters.7.md)).


# COMMANDS NOT IMPLEMENTED

`search_global`
: usearch's exhaustive database search, as opposed to the heuristic
  `usearch_global`. vsearch has no separate command for it because the
  heuristics can be switched off:
  `--usearch_global --maxaccepts 0 --maxrejects 0 --minwordmatches 0`
  compares every query to every target. The first two remove the early
  stop, the third removes the word pre-filter that selects candidates,
  and the alignment itself is always a full dynamic programming
  alignment. Expect the runtime that implies.

`cluster_otus`
: Not implemented, and there is no exact equivalent. It combines greedy
  clustering with chimera filtering (the UPARSE-OTU algorithm), so the
  closest vsearch route is a clustering command followed by a
  chimera-detection command, run separately.

`search_pcr`, `search_oligodb`
: Not implemented. Extracting the region between two primers is not
  something vsearch does; `--cut` cuts at a restriction pattern, which
  is a different operation (see
  [`vsearch-cut(1)`](../commands/vsearch-cut.1.md)).

`otutab_rare`, `otutab_norm`
: Not implemented. Rarefying or normalising an OTU table to a common
  number of reads is left to downstream tools. `--fastx_subsample`
  subsamples *sequences*, not table columns (see
  [`vsearch-fastx_subsample(1)`](../commands/vsearch-fastx_subsample.1.md)).

`otutab_octave`
: Not implemented.

`sintax_summary`
: Not implemented. `--sintax` writes per-query classifications; summing
  them per rank is left to downstream tools (see
  [`vsearch-sintax(1)`](../commands/vsearch-sintax.1.md)).


# OPTIONS SPELLED DIFFERENTLY

Same quantity, different name, and in one case a different direction.

| usearch                     | vsearch                            |
|-----------------------------|------------------------------------|
| `--fastq_pctid` *x*         | `--fastq_maxdiffpct` 100 - *x*     |
| `--nowordcountreject`       | `--minwordmatches 0`               |
| `--queryalnfract`           | `--query_cov`                      |
| `--targetalnfract`          | `--target_cov`                     |

`--fastq_pctid` bounds the identity of the overlap from below, where
`--fastq_maxdiffpct` bounds the percentage of mismatches from above, so
the two arguments are complements. Older usearch versions had
`--fastq_maxdiffpct` itself, under that name.


# OPTIONS ACCEPTED AND IGNORED

These options are recognised so that an existing command line still
runs, but they have no effect. `--band`, `--fulldp`, `--hspw`,
`--minhsp`, `--pattern`, `--slots` and `--xdrop_nw` all describe
usearch's alignment or seeding stage, which vsearch does not share:
vsearch always performs a full dynamic programming alignment, so there
is no band to widen and no heuristic to relax.

`--cons_truncate` is also accepted and ignored, with a warning.

`--strand plus` is accepted by `--uchime_ref` for compatibility, and is
the only value it takes; `--strand both` is a fatal error there.


# OPTIONS NOT RECOGNISED

Some options of usearch 5, 6 and 7, used by pipelines of that era such
as QIIME 1, are rejected outright rather than ignored:

`--global`
: vsearch always aligns globally, so the option would assert the only
  behaviour there is. Remove it.

`--evalue`
: vsearch does not compute E-values for nucleotide alignments. The
  `evalue` field of `--userfields` exists, and always reports -1 (see
  [`vsearch-userfields(7)`](./vsearch-userfields.7.md)).

`--query`
: Not an option in vsearch. The query file is the argument to the
  search command itself, as in `--usearch_global queries.fasta`. The
  name survives as the `query` field of `--userfields`.


# DELIBERATE DIFFERENCES

Corrections and extensions vsearch applies on purpose. They are
differences in output or in accepted input, not accidents.

- With `--blast6out` and `--output_no_hits`, usearch reports 13 fields
  for a query with no match, where the format has 12. vsearch reports
  the 12 the format calls for.
- With `--output_no_hits`, usearch lists queries without a match in the
  `--blast6out` file but not in the alignment output. vsearch lists them
  in both, as `No hits`.
- The `raw` field of `--userfields` is not informative in usearch.
  vsearch reports the alignment score.
- The fields `qlo`, `qhi`, `tlo` and `thi` have counterparts `qilo`,
  `qihi`, `tilo` and `tihi` reporting alignment coordinates that ignore
  terminal gaps.
- `--iddef`, and with it the alternative pairwise identity definitions
  that usearch removed, is available in vsearch.
- `--topn` extends to the sorting commands.
- `--sizein` extends to dereplication and clustering.
- T and U count as the same nucleotide during dereplication, so an RNA
  and a DNA spelling of one sequence collapse into a single entry.
- Sorting is stabilised, using sequence abundances or labels as
  secondary and tertiary keys, so equal keys do not come out in an
  arbitrary order.
- Low-complexity regions are masked with the DUST algorithm by default,
  and masking behaviour is more consistent (see
  [`vsearch-fastx_mask(1)`](../commands/vsearch-fastx_mask.1.md)).
- vsearch adds the command `--cluster_size`, which sorts sequences by
  decreasing abundance before clustering.

The `--tabbedout` report of `--fastq_mergepairs` follows the file usearch
writes for that option, with the same one-line-per-pair shape and the same
`result=merged` or `result=notmerged` ending each line. It differs in five
ways, all deliberate (see
[`vsearch-fastq_mergepairs(1)`](../commands/vsearch-fastq_mergepairs.1.md)):

- vsearch reports the percentage of *differences* in the overlap
  (`diffpct=`), where usearch reports the percentage of *identity*
  (`pctid=`). The two are complements, for the same reason `--fastq_pctid`
  and `--fastq_maxdiffpct` are (see OPTIONS SPELLED DIFFERENTLY above); the
  field is named after the option that bounds it.
- The three components of `aln=` are never negative. A pair whose reads run
  past each other is reported with a separate `stagger=` field giving the
  bases trimmed, where usearch encodes the same information as a negative
  first or third component of `aln=` -- which cannot be read by splitting the
  value on `-`.
- `trunc=` reports the read lengths that *remain* after `--fastq_truncqual`,
  where usearch's `tailf=` and `tailr=` report the number of bases *removed*.
- `relabel=` gives the label the merged sequence was written with, including
  any `;size=`, `;ee=` or `;length=` annotation requested. usearch has no
  such annotations, so its `relabel=` is always a bare label.
- Lines appear in input order whatever `--threads` is set to. usearch emits
  them in the order its threads finish, unless `-threads 1` is given.

usearch also reports a `minq=` field, for its `-fastq_minqual`. vsearch has
`--fastq_minqual`, but `--fastq_mergepairs` does not accept it (only
`--fastq_filter` and `--fastx_filter` do), so no such field is written.

vsearch also provides commands and options that usearch does not; they
are listed in [`vsearch(1)`](../index.1.md) and described one page per
command.


# SEE ALSO

[`vsearch(1)`](../index.1.md),
[`vsearch-nucleotides(7)`](./vsearch-nucleotides.7.md),
[`vsearch-pairwise_alignment_parameters(7)`](./vsearch-pairwise_alignment_parameters.7.md),
[`vsearch-userfields(7)`](./vsearch-userfields.7.md),
[`vsearch-usearch_global(1)`](../commands/vsearch-usearch_global.1.md)


#(../commands/fragments/footer.md)
