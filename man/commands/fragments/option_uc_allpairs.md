`--uc` *filename*
: Write pairwise alignment results to *filename* in a tab-separated
  uclust-like format with 10 columns. One line is written per alignment
  (record type always `H`). Columns are:

    1. record type, always `H`;
    2. ordinal number of the target sequence (zero-based);
    3. length of the query sequence;
    4. percentage of identity with the target;
    5. match orientation, always `+`;
    6. not used, always `0`;
    7. not used, always `0`;
    8. CIGAR alignment string (M, D, I; `=` if query and target are identical
       ignoring terminal gaps); see
       [`vsearch-cigar(5)`](../formats/vsearch-cigar.5.md);
    9. query label;
    10. target label.
