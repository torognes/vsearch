`--fastq_qmaxout` *integer*
: Specify the maximum quality score used when writing fastq files. The
  default is the highest score the output offset can represent (93 with
  `--fastq_asciiout` 33, 62 with offset 64), so scores read from the
  input are written back unchanged. Before version 3.0 the default was
  41, the usual maximum for Sanger/Illumina 1.8+ files, which silently
  reduced any higher score. Older formats may use a maximum quality
  score of 40. Two commands are exceptions and keep the old default of
  41, because they *generate* the score they clamp instead of passing
  one through: `--fasta2fastq`, which has no input quality and uses this
  option as the value to write, and `--fastq_mergepairs`, which caps the
  computed posterior quality of a merged base.

  For `--fastq_mergepairs` the offset in question is `--fastq_ascii`,
  not `--fastq_asciiout`: that command writes fastq but does not accept
  `--fastq_asciiout`, so a merged quality symbol carries the same offset
  the input was read with. The sum rule is stated against `--fastq_ascii`
  there, and against `--fastq_asciiout` everywhere else.
