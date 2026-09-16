`--uchimealns` *filename*
: Write three-way global alignments (parentA, parentB, chimera) to
  *filename* in a human-readable format. All sequences are converted
  to upper case before alignment. Lower case letters indicate
  disagreement in the alignment. Use `--alignwidth` to modify the
  alignment width. The alignment spans the parents over their full
  length, so a parent region reaching beyond the query is shown as a
  run of gaps on the query line rather than trimmed away.
