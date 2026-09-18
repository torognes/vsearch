`--denoise_omega_a` *real*
: Set the significance threshold for creating a new partition. A
  sequence becomes the center of a new partition when the abundance
  p-value of its reads, multiplied by the number of unique sequences
  (Bonferroni correction), is below this value. Lower values are more
  conservative: fewer, more abundant sequences are inferred. Accepted
  values range from 0.0 to 1.0 (excluded). The default is 1e-40, as
  in DADA2 (`OMEGA_A`).
