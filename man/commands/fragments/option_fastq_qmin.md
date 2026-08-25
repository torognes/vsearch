`--fastq_qmin` *integer*
: Specify the minimal quality score accepted when reading fastq
  sequences. Stop with an error message if a quality score lower than
  the specified value is read. The offset (see `--fastq_ascii`) plus
  the minimal score must be at least 33, the first printable ASCII
  character: scores down to 0 with offset 33, down to -31 with offset
  64. The value may therefore be negative, which is what the negative
  scores of older formats require, but note that those formats are not
  supported (see [`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md)).
  The default is 0, which is usual for recent Sanger/Illumina 1.8+
  files.
