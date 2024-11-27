`--fastq_asciiout` 33|64
: Specify the *offset* used as the basis for the fastq quality score
  when writing fastq output files. For example, an offset of 33 means
  that a quality value of 41 is represented by the 74th ASCII symbol
  (33 + 41 = 74), which is 'J'. See `ascii(7)` for a view of the ASCII
  character set. The offset value is either 33 or 64, default is 33.
