`--msaout` *filename*
: Write a multiple sequence alignment and a consensus sequence for each
  cluster to *filename*, in fasta format. vsearch computes center-star
  multiple sequence alignments using a fast method whose accuracy can
  decrease at low pairwise identity thresholds. Each cluster is
  written as the aligned centroid (its label prefixed with `*`),
  the aligned members, and a final `>consensus` record. In the
  consensus, each column within the centroid's span shows its most
  frequent nucleotide, or `-` when gaps outnumber every nucleotide;
  columns outside the centroid's span (terminal extensions contributed
  by longer members) are shown as `+`, whatever their content. If
  `--sizein` is specified, sequence abundances are taken into account
  when computing the consensus.
