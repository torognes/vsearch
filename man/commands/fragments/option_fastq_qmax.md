`--fastq_qmax` *integer*
: Specify the maximal quality score accepted when reading fastq
  sequences. Stop with an error message if a quality score higher than
  the specified value is read. The offset (see `--fastq_ascii`) plus
  the maximal score may not exceed 126, the last printable ASCII
  character: scores up to 93 with offset 33, up to 62 with offset
  64. The default is the highest score the offset can represent (93
  with offset 33, 62 with offset 64), so no quality score is rejected
  unless this option is lowered. Before version 3.0 the default was
  41, the usual maximum for Sanger/Illumina 1.8+ files, which rejected
  PacBio HiFi and nanopore files outright.
