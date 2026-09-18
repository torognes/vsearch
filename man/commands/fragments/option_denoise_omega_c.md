`--denoise_omega_c` *real*
: Set the threshold below which a read is not corrected. A unique
  sequence whose abundance p-value (unconditional) within its final
  partition is below this value is too abundant to be an error of the
  center, yet was not significant enough to found a partition: its
  reads are left out of `--fastqout` and written to
  `--fastqout_discarded`. With 0, all reads are corrected. Accepted
  values range from 0.0 to 1.0 (excluded). The default is 1e-40, as
  in DADA2 (`OMEGA_C`).
