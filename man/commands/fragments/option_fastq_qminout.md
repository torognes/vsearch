`--fastq_qminout` *integer*
: Specify the minimum quality score used when writing fastq files. The
  offset (see `--fastq_asciiout`) plus the minimum score must be at
  least 33, so the value may be negative when the output offset is 64.
  The default is 0, which is usual for recent Sanger/Illumina 1.8+
  files. Older formats may use scores between -5 and 2.

  For `--fastq_mergepairs` the offset in question is `--fastq_ascii`,
  not `--fastq_asciiout`: that command writes fastq but does not accept
  `--fastq_asciiout`, so a merged quality symbol carries the same offset
  the input was read with, and the bound follows it. With
  `--fastq_ascii` 64 the minimum score may therefore go down to -31.
