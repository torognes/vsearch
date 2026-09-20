`--tabbedout` *filename*
: Write to *filename* a report of how each pair of reads fared, one line
  per input pair, whether or not it merged. Pairs appear in input order,
  whatever the value of `--threads`. Each line is a sequence of
  tab-separated tokens: the label of the forward read, then a token for
  each value the merging pipeline computed, then the reason the pair was
  rejected if it was, then `result=merged` or `result=notmerged`. A value
  token is present only when the stage producing it was reached, so a pair
  rejected early yields a short line. For example:

    ~~~
    s1	len=20-20	aln=0-20-0	diffs=1	diffpct=5.0	mergelen=20	ee=0.5027	relabel=s1	result=merged
    s2	len=20-20	aln=0-20-0	diffs=1	diffpct=5.0	mergelen=20	toomanydiffs	result=notmerged
    s3	len=1-1	nokmers	result=notmerged
    ~~~

  The value tokens are:

    - `len=`*forward*`-`*reverse*: the lengths of the two input reads,
      before any truncation. Always present.
    - `trunc=`*forward*`-`*reverse*: the lengths that remain after
      `--fastq_truncqual` has trimmed low-quality 3' tails. Present only
      when that trimming shortened at least one of the two reads.
    - `aln=`*left*`-`*overlap*`-`*right*: the alignment, as the number of
      bases the forward read contributes before the overlap, the length of
      the overlap, and the number of bases the reverse read contributes
      after it. The three values are never negative and always sum to
      `mergelen=`. Present only when an alignment was found.
    - `stagger=`*forward*`-`*reverse*: the number of 3' bases of each read
      that run past the other read's 5' end. Such bases are trimmed when
      the pair is merged (see `--fastq_allowmergestagger`). Present only
      when at least one of the two is non-zero.
    - `diffs=`*integer*: the number of mismatches in the overlap.
    - `diffpct=`*float*: those mismatches as a percentage of the overlap
      length, the quantity `--fastq_maxdiffpct` bounds.
    - `mergelen=`*integer*: the length of the merged sequence the
      alignment implies.
    - `ee=`*float*: the number of expected errors in the merged sequence.
      Present only once the merged sequence has been assembled. See
      [`vsearch-expected_error(7)`](../misc/vsearch-expected_error.7.md).
    - `relabel=`*label*: the label the merged sequence was written with,
      including any annotation added by `--relabel`, `--sizeout`,
      `--fastq_eeout` or `--lengthout`. Present only for merged pairs, as
      nothing was written for the others. The first field of the line
      always remains the input label.

  When a pair is not merged, exactly one of the following tokens names the
  reason, immediately before `result=notmerged`:

    - `tooshort`: a read is shorter than `--fastq_minlen`;
    - `toolong`: a read is longer than `--fastq_maxlen`;
    - `toomanyns`: a read has more N's than `--fastq_maxns` allows;
    - `nokmers`: too few k-mers were found on the same diagonal, so no
      alignment could be attempted;
    - `multiplealns`: several alignments of comparable quality were found,
      suggesting a repeated region;
    - `nostagger`: the pair is staggered and `--fastq_allowmergestagger`
      was not specified;
    - `toomanydiffs`: more mismatches than `--fastq_maxdiffs` allows;
    - `toomanydiffpct`: a higher percentage of mismatches than
      `--fastq_maxdiffpct` allows;
    - `lowscore`: the alignment score is too low, or its drop too steep,
      which usually points to clustered mismatches or an indel in the
      overlap region;
    - `alntooshort`: the overlap is shorter than `--fastq_minovlen`;
    - `mergetooshort`: the merged sequence is shorter than
      `--fastq_minmergelen`;
    - `mergetoolong`: the merged sequence is longer than
      `--fastq_maxmergelen`;
    - `toohighee`: the merged sequence has more expected errors than
      `--fastq_maxee` allows.

  This report is modelled on the file usearch writes for the same option,
  but it is not identical to it; the differences are deliberate and are
  listed in [`vsearch-usearch(7)`](../misc/vsearch-usearch.7.md).
