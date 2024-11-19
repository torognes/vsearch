`--xee`
: Strip expected error (ee) annotations from fastq or fasta headers
  when writing output files (find and remove the pattern
  `;ee=float[;]` from sequence headers). Expected error annotations
  are added by the synonymous options `--fastq_eeout` and `--eeout`
  (see for example command `--fastx_filter`, described in
  [`vsearch-fastx_filter(1)`](./commands/vsearch-fastx_filter.1.md)).
