`--fastq_qmax` *positive integer*
: Specify the maximal quality score accepted when reading fastq
  sequences. Reads with higher quality scores are discarded. Accepted
  values range from 0 to 93 if the offset is 33 (see `--fastq_ascii`),
  or range from 0 to 62 if the offset is 64. The default is 41, which
  is usual for recent Sanger/Illumina 1.8+ files.
