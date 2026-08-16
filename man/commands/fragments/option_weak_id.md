`--weak_id` *real*
: Report hits with a pairwise identity of at least *real*, without
  stopping the search. Unlike `--id`, weak hits do not count toward
  `--maxaccepts` but do count toward `--maxrejects`. Values larger
  than the value specified with `--id` are silently reduced to it.
  With `--cluster_unoise`, which does not use `--id`, *real* is the
  identity floor of the denoising step and defaults to 0.90.
