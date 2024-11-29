`--fastq_qmin` *positive integer*
: Specify the minimal quality score accepted when reading fastq
  sequences. Reads with lower quality scores are discarded. Accepted
  values range from 0 to 93 if the offset is 33 (see `--fastq_ascii`),
  or range from 0 to 62 if the offset is 64. The default is 0, which
  is usual for recent Sanger/Illumina 1.8+ files. Older formats may
  use scores between -5 and 2.
