`--fastq_qmax` *positive integer*
: Specify the maximal quality score accepted when reading fastq
  sequences. Stop with an error message if a quality score higher than
  the specified value is read. The offset (see `--fastq_ascii`) plus
  the maximal score may not exceed 126, the last printable ASCII
  character: scores up to 93 with offset 33, up to 62 with offset
  64. The default is 41, which is usual for recent Sanger/Illumina
  1.8+ files.
