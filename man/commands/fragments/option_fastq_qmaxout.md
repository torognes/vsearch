`--fastq_qmaxout` *integer*
: Specify the maximum quality score used when writing fastq files. The
  default is the highest score the output offset can represent (93 with
  `--fastq_asciiout` 33, 62 with offset 64), so scores read from the
  input are written back unchanged. Before version 3.0 the default was
  41, the usual maximum for Sanger/Illumina 1.8+ files, which silently
  reduced any higher score. Older formats may use a maximum quality
  score of 40. The sole exception is `--fasta2fastq`, which has no
  input quality and uses this option as the value to write rather than
  as a ceiling: there the default remains 41.
