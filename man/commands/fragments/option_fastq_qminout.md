`--fastq_qminout` *integer*
: Specify the minimum quality score used when writing fastq files. The
  offset (see `--fastq_asciiout`) plus the minimum score must be at
  least 33, so the value may be negative when the output offset is 64.
  The default is 0, which is usual for recent Sanger/Illumina 1.8+
  files. Older formats may use scores between -5 and 2.

  For `--fastq_mergepairs` the offset in question is `--fastq_ascii`,
  not `--fastq_asciiout`: that command writes fastq but does not accept
  `--fastq_asciiout`, so a merged quality symbol carries the same offset
  the input was read with, and the bound follows it.

  A negative minimum matters where the score passes through vsearch and
  may itself be negative: `--fastq_convert` on an offset-64 file read
  with a lowered `--fastq_qmin` writes the score back unchanged with
  `--fastq_qminout` -5, where the default 0 would raise it. It has no
  effect where vsearch *computes* the score it clamps, because such a
  score is derived from an error probability and is therefore never
  negative: `--fastq_mergepairs`, and `--fastq_convert --fastq_solexa`,
  which converts to the Phred scale before any output clamp applies.
