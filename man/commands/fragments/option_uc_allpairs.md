`--uc` *filename*
: Write pairwise alignment results to *filename* in a tab-separated
  uclust-like format with 10 columns. One line is written per
  alignment (record type `H`); a query with no accepted alignment
  produces a no-hit record instead (record type `N`) --- the last
  sequence of the file always does, as it has no following targets to
  be aligned with. Columns are:

    1. record type: `H` (alignment) or `N` (no alignment);
    2. ordinal number of the target sequence (zero-based; `*` for N);
    3. length of the query sequence (`*` for N);
    4. percentage of identity with the target (`*` for N);
    5. match orientation `+` (`.` for N);
    6. not used, always `0` (`*` for N);
    7. not used, always `0` (`*` for N);
    8. CIGAR alignment string (M, D, I; `=` if the query and target
       sequences are strictly identical; `*` for N); see
       [`vsearch-cigar(5)`](../formats/vsearch-cigar.5.md);
    9. query label;
    10. target label (`*` for N).
