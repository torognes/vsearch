`--uc` *filename*
: Write search results to *filename* in a tab-separated uclust-like format
  with 10 columns. Each query produces a hit (`H`) or no-hit (`N`) record.
  Use `--uc_allhits` to report all hits per query instead of just the top
  hit. Columns are:

    1. record type: `H` (hit) or `N` (no hit);
    2. ordinal number of the target sequence (zero-based; `*` for N);
    3. length of the query sequence (`*` for N);
    4. percentage of identity with the target (`*` for N);
    5. match orientation `+` or `-` (`.` for N);
    6. not used; `0` (H) or `*` (N);
    7. not used; `0` (H) or `*` (N);
    8. CIGAR alignment string (`*` for N);
    9. query label;
    10. target label (`*` for N).
