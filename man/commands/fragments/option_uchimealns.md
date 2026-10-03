`--uchimealns` *filename*
: Write three-way global alignments (parentA, parentB, chimera) to
  *filename* in a human-readable format. All sequences are converted
  to upper case before alignment. Lower case letters indicate
  disagreement in the alignment. Use `--alignwidth` to modify the
  alignment width. The alignment is trimmed to the query: parent
  regions reaching beyond either end of the query are not shown, and
  the positions on the parent lines count from the start of each
  parent (so a parent line can start at a position greater than 1).
  These regions are not used to detect chimeras.
