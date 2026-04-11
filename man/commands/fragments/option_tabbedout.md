`--tabbedout` *filename*
: Write dereplication details to *filename* as a tab-separated file
  with one row per input sequence and 6 columns:

    1. original sequence label;
    2. output label (label of the first sequence in the cluster,
       possibly relabelled);
    3. cluster number (zero-based);
    4. sequence number within the cluster (zero-based);
    5. cluster size;
    6. original label of the first sequence in the cluster (before
       any relabelling).
