% vsearch-usearch(7) version 2.33.0 | vsearch manual
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
implement, commands and options it spells differently, options it
accepts and ignores, options it rejects, and behaviours it changed on
purpose.

Two differences are broad enough to state first. vsearch works with
nucleotide sequences only; amino acid sequences are not supported (see
[`vsearch-nucleotides(7)`](./vsearch-nucleotides.7.md)). And vsearch
always aligns globally, with full dynamic programming, so the options
by which usearch tunes or bounds its alignment stage have nothing to
act on (see
[`vsearch-pairwise_alignment_parameters(7)`](./vsearch-pairwise_alignment_parameters.7.md)).


# COMMANDS NOT IMPLEMENTED

`cluster_otus`
: Not implemented, and there is no exact equivalent. It combines greedy
  clustering with chimera filtering (the UPARSE-OTU algorithm), so the
  closest vsearch route is a clustering command followed by a
  chimera-detection command, run separately. usearch's own manual
  considers 97% OTU clustering obsolete for most purposes, and
  recommends denoising instead --- a route vsearch does provide (see
  `unoise3` in the next section).

`search_pcr`, `search_pcr2`, `search_oligodb`
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


# COMMANDS SPELLED DIFFERENTLY

The operation is available, under another name, or as an option of
another command.

`unoise3`
: `--cluster_unoise`, then `--uchime3_denovo`. usearch denoises and
  removes chimeras in a single command; vsearch splits the two, so the
  chimera step is a separate run (see
  [`vsearch-cluster_unoise(1)`](../commands/vsearch-cluster_unoise.1.md)).

`otutab`
: `--usearch_global` with `--otutabout`. The OTU sequences that usearch
  passes with `-otus` or `-zotus` are the `--db` here, and the identity
  threshold that usearch fixes at 0.97 by default is the usual `--id`
  (see
  [`vsearch-usearch_global(1)`](../commands/vsearch-usearch_global.1.md)).
  `--search_exact` also writes `--otutabout`, when only exact matches
  should count.

`fastx_truncate`
: `--fastx_filter`, whose `--fastq_stripleft`, `--fastq_stripright` and
  `--fastq_trunclen` perform the same trimming; `--fastq_trunclen_keep`
  keeps the sequences that are shorter than the requested length,
  where usearch discards them. Padding (`-padlen`, `-padq`) has no
  equivalent (see
  [`vsearch-fastx_filter(1)`](../commands/vsearch-fastx_filter.1.md)).

`fastx_orient`
: `--orient`. usearch 11 accepts both spellings, `orient` and
  `fastx_orient`; usearch 12 keeps only the second (see
  [`vsearch-orient(1)`](../commands/vsearch-orient.1.md)).

`cluster_fast` with `-sort size`
: `--cluster_size`. vsearch names the sort order rather than passing
  it: `--cluster_fast` sorts by decreasing length, `--cluster_size` by
  decreasing abundance, and `--cluster_smallmem` does not sort at all,
  expecting its input already sorted by decreasing length (see
  [`vsearch-cluster_size(1)`](../commands/vsearch-cluster_size.1.md)).


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

- `--search_global` ignores `--maxaccepts` and `--maxrejects`, which
  usearch honours. An exhaustive search has no ranked candidate order
  for an early stop to cut short, so honouring them would make the
  reported hits depend on how the database happens to be sorted, which
  is the one property the command exists to remove. They are accepted
  and have no effect, so a command line written for `--usearch_global`
  can be reused unchanged; `--maxhits` limits how many hits are
  reported. `--search_global` also defaults `--minseqlength` to 1 rather
  than 32, nothing in it requiring a sequence to hold a whole word (see
  [`vsearch-search_global(1)`](../commands/vsearch-search_global.1.md)).
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
- This includes `--sintax`, which DUST-masks a fasta reference database
  by default (`--dbmask dust`). usearch's `-sintax` masks nothing unless
  `-dbmask` is given, so low-complexity reference regions can take part
  in its classifications and not in vsearch's; `--dbmask none` restores
  usearch's default. The two also treat lowercase letters differently:
  vsearch's *dust* ignores the case of the input and masks only
  low-complexity regions, whereas usearch's `-dbmask dust` keeps any
  lowercase region masked as well, so an all-lowercase reference
  database leaves every query unclassified in usearch
  (see [`vsearch-sintax(1)`](../commands/vsearch-sintax.1.md)).
- vsearch adds the command `--cluster_size`, which sorts sequences by
  decreasing abundance before clustering.
- `--relabel @`, which derives a sample identifier from the input file
  name, is accepted by every command that accepts `--relabel`. usearch
  honours `-relabel @` for `fastq_filter` and `fastq_mergepairs` only,
  and takes it as a literal `@` prefix everywhere else, so a command
  such as `--fastx_uniques` writes '>sampleA.1' where usearch writes
  '>@1'. vsearch also refuses `--relabel @` on standard input, which has
  no file name to derive from, and warns when the name yields an empty
  identifier; usearch cannot read standard input at all. The derivation
  itself follows usearch's implementation rather than its documentation
  (see `--relabel` in any command's manual page).
- The output-label column of `--fastx_uniques --tabbedout` follows
  `--relabel_md5`, `--relabel_self` and `--relabel_sha1` as well as
  `--relabel`. usearch has none of those three options, so the column is
  byte-identical to its own whenever `-relabel` is what was given.

The `--tabbedout` report of `--fastq_mergepairs` follows the file usearch
writes for that option, with the same one-line-per-pair shape and the same
`result=merged` or `result=notmerged` ending each line. It differs in the
following ways, all deliberate (see
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
- The first field is the whole header of the forward read, description
  included, which is the label the merged-read outputs of the command carry.
  usearch truncates its first field at the first blank.

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
