`--fastq_solexa`
: Convert quality scores from the Solexa/Illumina 1.0 scale to the
  Phred scale while re-encoding. A Solexa score is *Q = -10 log10(p /
  (1 - p))*, not the Phred *Q = -10 log10(p)*, so rebasing the ASCII
  offset alone would leave the scores wrong (see
  [`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md)). The conversion
  applied is *Q_phred = 10 log10(10^(Q_solexa / 10) + 1)*, rounded to
  the nearest integer; it is the identity from Solexa 10 upward, so
  only the fifteen lowest scores change. Requires `--fastq_ascii 64`,
  the offset the Solexa encoding uses, and implies `--fastq_qmin -5`,
  the floor of the Solexa scale, unless `--fastq_qmin` is given
  explicitly. The conversion is one-way and lossy: six Solexa pairs
  collapse onto a single Phred score each, and there is no option to
  write Solexa scores. Only `--fastq_convert` accepts this option.
