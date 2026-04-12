`--profile` *filename*
: Write a sequence profile to *filename*, reporting the frequency of
  each nucleotide at each position in the multiple alignment for each
  cluster. A FASTA-like header line precedes the profile information
  for each cluster. The data is tab-separated with eight columns:
  position (0-based), consensus nucleotide, number of As, Cs, Gs, Ts
  or Us, gap symbols, and total number of ambiguous nucleotide symbols
  (B, D, H, K, M, N, R, S, Y, V or W). If `--sizein` is specified,
  sequence abundances are taken into account.
