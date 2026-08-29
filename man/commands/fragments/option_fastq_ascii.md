`--fastq_ascii` 33|64
: Specify the *offset* used as the basis for the fastq quality score
  when reading fastq files. For example, an offset of 33 means that a
  quality value of 41 is represented by the 74th ASCII symbol (33 + 41
  = 74), which is 'J'. See `ascii(7)` for a view of the ASCII
  character set. The offset value is either 33 or 64, default is 33.

  The offset matters even to a command that never decodes a quality
  score: it is what the reader compares the observed quality symbols
  against before warning that the file may use the other encoding, it
  sets the default of `--fastq_qmax` (the highest score the offset can
  represent), and the sum rules on `--fastq_qmin` and `--fastq_qmax` are
  stated in terms of it.
