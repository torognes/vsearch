`--mothur_shared_out` *filename*
: Write an OTU table to *filename* in the mothur 'shared'
  tab-separated plain text format. The first line starts with `label`,
  `group` and `numOtus`, followed by all OTU identifiers. Each
  subsequent line starts with `vsearch`, the sample identifier, the
  total number of OTUs, and the abundance of each OTU in that sample.
  Sample and OTU identifiers are extracted from FASTA headers. OTUs
  are represented by the cluster centroids. Abundance annotations
  (`;size=integer`) present in sequence headers are always used when
  filling the table, whether or not `--sizein` is given (unlike the
  `--uc` cluster summaries, which count each sequence as 1 without
  `--sizein`).
  OTU identifiers must be unique: two OTUs sharing the same identifier
  are reported in a single column, `numOtus` counts that column only
  once, and their abundances are summed. When clustering, a relabelling
  option (`--relabel`, `--relabel_self`, `--relabel_md5` or
  `--relabel_sha1`) guarantees unique identifiers; when searching, the
  database itself must have unique headers.
