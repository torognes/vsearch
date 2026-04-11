`--uc` *filename*
: Write dereplication results to *filename* in a tab-separated
  uclust-like format with 10 columns. Three record types are used per
  row: cluster seeds (`S`), hits (`H`), and cluster summaries (`C`).
  Columns are:

    1. record type (`S`, `H`, or `C`);
    2. cluster number (zero-based);
    3. sequence length (`S`, `H`), or cluster size (`C`);
    4. percent identity with centroid (`H`), or `*` (`S`, `C`);
    5. match orientation `+` or `-` (`H`), or `*` (`S`, `C`);
    6–8. unused, always `*` (`S`, `C`) or `0` (`H`);
    9. query label (`H`), or centroid label (`S`, `C`);
    10. centroid label (`H`), or `*` (`S`, `C`).
