`--dn` *strictly positive real number*
: Set the pseudo-count prior on the number of 'no' votes,
  corresponding to the parameter *n* in the chimera scoring function.
  Default value is 1.4. Increasing `--dn` reduces the likelihood of
  tagging a sequence as a chimera (fewer false positives, but also
  more false negatives).
