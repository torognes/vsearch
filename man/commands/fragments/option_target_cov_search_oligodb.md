`--target_cov` *real*
: Report an occurrence truncated by a read end only if at least
  *real* of the oligo lies inside the read. The default is 0.75. The
  oligo nucleotides hanging off the read end cost nothing (they are not
  counted in `diffs`), but the shorter the part left inside the read,
  the more often it matches by chance: at 0.5, random partial matches
  outnumber real truncated oligos (see DESCRIPTION). Use
  `--target_cov 1` to report whole occurrences only.

    Unlike the other commands, where the target coverage counts only
    the target nucleotides aligned with a query nucleotide, the
    fraction compared against here is the span of the oligo inside the
    read, deleted nucleotides included: a deletion is a difference,
    bounded by `--maxdiffs`, not a truncation. The `tcov` userfield
    keeps its usual definition, so it is below 100% for an occurrence
    with a deletion even when the whole oligo lies inside the read;
    `tlo` > 1 or `thi` < `tl` tells a truncated occurrence apart.
