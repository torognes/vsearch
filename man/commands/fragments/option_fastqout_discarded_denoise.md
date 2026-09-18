`--fastqout_discarded` *filename*
: Write the reads that were not corrected to *filename*, unchanged, in
  fastq format: reads with ambiguous nucleotides (N), reads of five
  nucleotides or less, and reads rejected by `--denoise_omega_c`.
