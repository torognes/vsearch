`--tabbedout` *filename*
: Write classification results to *filename* as a tab-separated file. Each
  query produces one row with the following columns:

    1. query label;
    2. predicted taxonomy with bootstrap confidence in parentheses after each
       rank (e.g., `d:Bacteria(1.00),p:Proteobacteria(0.95)`);
    3. strand (`+` or `-`);
    4. predicted taxonomy filtered by `--sintax_cutoff`, showing only ranks
       with bootstrap support at or above the threshold and omitting the
       values (e.g., `d:Bacteria,p:Proteobacteria`). Absent if `--sintax_cutoff`
       is not specified.

  This option is mandatory.
