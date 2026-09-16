`--otutabout` *filename*
: Write an OTU table to *filename* in a classic tab-separated plain
  text format. The first line starts with `#OTU ID` followed by sample
  identifiers. Each subsequent line starts with the OTU identifier
  followed by the abundances in each sample. Sample and OTU
  identifiers are extracted from FASTA headers (see `--sample`). OTUs
  are represented by the cluster centroids. A `taxonomy` column is
  appended if taxonomy information is available for any OTU. Abundance
  annotations (`;size=integer`) present in sequence headers are always
  used when filling the table, whether or not `--sizein` is given
  (unlike the `--uc` cluster summaries, which count each sequence as 1
  without `--sizein`).
  Rows and columns are both ordered lexicographically by identifier,
  whatever order the OTUs were created in or the samples first appeared
  in, and no option changes that. Sorting the table by decreasing
  abundance is a job for a downstream tool.

  OTU identifiers must be unique: two OTUs sharing the same identifier
  are reported on a single line, and their abundances are summed. When
  clustering, a relabelling option (`--relabel`, `--relabel_self`,
  `--relabel_md5` or `--relabel_sha1`) guarantees unique identifiers;
  when searching, the database itself must have unique headers.
