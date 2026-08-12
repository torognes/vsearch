`--strand` *plus*
: Check the *plus* strand only (default). Unlike full-length
  dereplication, `--derep_prefix` does not accept `--strand both` and
  rejects it with a fatal error: reverse-complementation does not
  preserve the prefix relation (a prefix on the minus strand is a
  suffix on the plus strand), so prefix groups would depend on the
  input order. Reads in mixed orientations should be oriented first
  (see [`vsearch-orient(1)`](./vsearch-orient.1.md)).
