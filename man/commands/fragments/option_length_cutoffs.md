`--length_cutoffs` *min*,*max*,*increment*
: Set the range of truncation lengths to use as rows in the output
  table. Requires three comma-separated integers: the shortest cutoff,
  the longest cutoff, and the increment between cutoffs. The longest
  cutoff may be `*` to use the length of the longest sequence in the
  input file. The default is `50,*,50`, producing rows at 50, 100,
  150, and so on up to the longest sequence.
