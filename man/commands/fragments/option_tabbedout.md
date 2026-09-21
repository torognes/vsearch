`--tabbedout` *filename*
: Write dereplication details to *filename* as a tab-separated file
  with one row per input sequence and 6 columns:

    1. original sequence label;
    2. output label: the label of the first sequence in the cluster,
       or the substitute requested by `--relabel`, `--relabel_md5`,
       `--relabel_self` or `--relabel_sha1`. This is the bare label:
       the annotations `--sizeout`, `--lengthout` and `--sample` add to
       the fasta header are not repeated here;
    3. cluster number (zero-based);
    4. sequence number within the cluster (zero-based);
    5. cluster size (the number of sequences merged into the
       cluster, whether or not `--sizein` is used);
    6. original label of the first sequence in the cluster (before
       any relabelling).
