# VSEARCH5D

This is a modified version of VSEARCH. See https://github.com/torognes/vsearch to know about VSEARCH.

## Difference from VSEARCH

The --idoffset argument was added. If you specify --idoffset 10, the identity will be calculated based on the following command.

100 * (the number of matches - idoffset) / (alignment length - idoffset)

This argument should be useful for clustering of non-overlapped paired-end reads as explained below.

Assuming you constructed concatenated sequences like the following,
> <- reverse_read - padding - forward_read ->
the clustered results will be OVERclustered because of padding sequence. In such cases, giving the length of padding sequence as idoffset will produce correctly clustered results. After the clustering, the representative or consensus sequences can be divided into forward and reverse parts using padding sequence.
