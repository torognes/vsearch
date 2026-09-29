`--maxgaps` *positive integer*
: Report an occurrence only if its alignment has at most *integer* gap
  openings: a run of consecutive gap columns counts as one, whatever
  its length (`--maxdiffs` bounds the number of gap columns). The
  default is unlimited. The alignment compared against `--maxgaps` is
  the one with the fewest differences, and among those, the one with
  the fewest gap openings; an occurrence is rejected when that
  alignment has too many gap openings, even if an alignment without
  gaps but with more substitutions would fit within `--maxdiffs`.

    `--maxgaps 0` is the exception: it switches to substitutions only.
    The search then looks for occurrences with up to `--maxdiffs`
    mismatches and no gap at all, as usearch does, and gapped
    alignments are not considered.
