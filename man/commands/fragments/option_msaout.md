`--msaout` *filename*
: Write a multiple sequence alignment and a consensus sequence for each
  cluster to *filename*, in fasta format. vsearch computes center-star
  multiple sequence alignments using a fast method whose accuracy can
  decrease at low pairwise identity thresholds. The consensus sequence
  is constructed by taking the majority symbol (nucleotide or gap) from
  each column of the alignment. Columns containing a majority of gaps
  are skipped, except for terminal gaps. If `--sizein` is specified,
  sequence abundances are taken into account when computing the
  consensus.
