`--allow_fewer`
: Treat the `--sample_size` value as an upper bound: subsample that
  number of reads, or output all of them when fewer are available,
  instead of reporting a fatal error. The reduced output is expected
  and is not signalled in any way. This option is ignored when the
  subset size is given as a percentage with `--sample_pct`, which can
  never exceed the input; a warning is issued in that case.
