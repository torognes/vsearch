`--uc` *filename*
: Write results to *filename* in a tab-separated uclust-like format
  with 10 columns. Three record types are used per row: cluster seeds
  (`S`), hits (`H`), and cluster summaries (`C`). Columns are:

    1. record type (`S`, `H`, or `C`);
    2. cluster number (zero-based);
    3. centroid length (`S`), query length (`H`), or cluster size (`C`);
    4. percent identity with centroid (`H`), or `*` (`S`, `C`);
    5. match orientation `+` or `-` (`H`), or `*` (`S`, `C`);
    6. not used; always `0` (`H`) or `*` (`S`, `C`);
    7. not used; always `0` (`H`) or `*` (`S`, `C`);
    8. CIGAR alignment string (`H`), or `*` (`S`, `C`);
    9. query label (`H`), or centroid label (`S`, `C`);
    10. centroid label (`H`), or `*` (`S`, `C`).
