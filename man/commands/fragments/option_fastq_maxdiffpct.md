`--fastq_maxdiffpct` *real*
: Set the maximum percentage of mismatches allowed in the overlap
  region (0.0 to 100.0). Additional heuristics in the merging
  algorithm may still discard pairs with a high mismatch rate. The
  default is 100.0. usearch bounds the same quantity from the other
  side, as a minimum identity: `--fastq_maxdiffpct` *x* is equivalent
  to usearch's `--fastq_pctid` 100 - *x*.
