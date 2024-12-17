% vsearch-fastq_stats(1) version 2.30.0 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./fragments/date.md)

# NAME

vsearch \-\-fastq_stats --- analyze fastq sequences and output
detailed statistics


# SYNOPSIS

| **vsearch** **\-\-fastq_stats** _fastqfile_ \-\-log _outputfile_ \[_options_]


# DESCRIPTION

The vsearch command `--fastq_stats` analyzes fastq sequences and
outputs detailed statistics, including projections of the impact
different filtering thresholds would have. Results are written to
`--log` *outputfile*. The quality encoding of the input fastq and the
range of accepted quality values can be specified with
`--fastq_ascii`, `--fastq_qmin` and `--fastq_qmax`.

The five different statistics tables reported by `--fastq_stats` are
described below, with an example for each.

## Read length distribution

Observed read lengths are sorted in decreasing order. The largest
length value is marked with '>='.

1.  L: read length.
2.  N: number of reads.
3.  Pct: fraction of reads with this length.
4.  AccPct: fraction of reads with this length or longer.

```text
Read length distribution
      L           N      Pct   AccPct
-------  ----------  -------  -------
>=  251       19795    69.2%    69.2%
    250        6847    23.9%    93.2%
    249        1626     5.7%    98.9%
    248         260     0.9%    99.8%
    247          24     0.1%    99.8%
[...]
     75           5     0.0%   100.0%
     74           2     0.0%   100.0%
     66           2     0.0%   100.0%
     65           4     0.0%   100.0%
     64           1     0.0%   100.0%
```


## Quality score distribution

Observed quality values are sorted in decreasing order.

1.  ASCII: character encoding the quality score.
2.  Q: Phred quality score.
3.  Pe: probability of error associated with the quality score.
4.  N: number of bases with this quality score.
5.  Pct: fraction of bases with this quality score.
6.  AccPct: fraction of bases with this quality score or higher.

```text
Q score distribution
ASCII    Q       Pe           N      Pct   AccPct
-----  ---  -------  ----------  -------  -------
    I   40  0.00010           2     0.0%     0.0%
    H   39  0.00013     1009552    14.1%    14.1%
    G   38  0.00016     1236861    17.3%    31.4%
    F   37  0.00020     1231882    17.2%    48.6%
    E   36  0.00025      330094     4.6%    53.2%
    D   35  0.00032      143134     2.0%    55.2%
    C   34  0.00040      250511     3.5%    58.7%
    B   33  0.00050      451696     6.3%    65.0%
    A   32  0.00063      273153     3.8%    68.8%
    @   31  0.00079      171356     2.4%    71.2%
    ?   30  0.00100      131957     1.8%    73.0%
    >   29  0.00126       59308     0.8%    73.9%
    =   28  0.00158       25174     0.4%    74.2%
    <   27  0.00200       44282     0.6%    74.8%
    ;   26  0.00251      204896     2.9%    77.7%
    :   25  0.00316       49600     0.7%    78.4%
    9   24  0.00398      251427     3.5%    81.9%
    6   21  0.00794          32     0.0%    81.9%
    5   20  0.01000       11170     0.2%    82.0%
    4   19  0.01259        9188     0.1%    82.2%
    3   18  0.01585       34528     0.5%    82.7%
    2   17  0.01995       25452     0.4%    83.0%
    1   16  0.02512      100936     1.4%    84.4%
    0   15  0.03162       96675     1.3%    85.8%
    /   14  0.03981      387620     5.4%    91.2%
    .   13  0.05012      162798     2.3%    93.5%
    -   12  0.06310      468864     6.5%   100.0%
```


## Length vs. quality distribution

Positions in reads are sorted in increasing order.

1.  L: position in reads (starting from position 2).
2.  PctRecs: fraction of reads with at least this length.
3.  AvgQ: average quality at this position.
4.  P(AvgQ): error probability corresponding to AvgQ.
5.  AvgP: average error probability at this position.
6.  AvgEE: average expected error over all reads up to this position.
7.  Rate: growth rate of AvgEE between the current and previous position.
8.  RatePct: Rate (as explained above) expressed as a percentage.

```text
    L  PctRecs  AvgQ  P(AvgQ)      AvgP  AvgEE       Rate   RatePct
-----  -------  ----  -------  --------  -----  ---------  --------
    2   100.0%  28.2  0.00153  0.005655   0.01   0.005053    0.505%
    3   100.0%  30.1  0.00097  0.002893   0.01   0.004333    0.433%
    4   100.0%  27.7  0.00171  0.005935   0.02   0.004734    0.473%
    5   100.0%  28.8  0.00131  0.004416   0.02   0.004670    0.467%
    6   100.0%  29.5  0.00112  0.005633   0.03   0.004830    0.483%
[...]
  247    99.8%  24.4  0.00360  0.019798   2.14   0.008668    0.867%
  248    99.8%  24.0  0.00394  0.020644   2.16   0.008711    0.871%
  249    98.9%  23.9  0.00407  0.020943   2.17   0.008716    0.872%
  250    93.2%  23.8  0.00417  0.021239   2.22   0.008877    0.888%
  251    69.2%  22.3  0.00589  0.023948   2.27   0.009048    0.905%
```


## Effect of expected error and length filtering

The first column indicates read lengths (L), sorted in decreasing
order, and starting with the first L with a cummulated expected error
smaller or equal to 1.0. The next four columns indicate the number of
reads that would be retained by the `--fastq_filter` command if the
reads were truncated at length L (option `--fastq_trunclen L`) and
filtered to have a maximum expected error of 1.0, 0.5, 0.25 or 0.1
(with the option `--fastq_maxee float`). The last four columns
indicate the fraction of reads that would be retained by the
`--fastq_filter` command using the same length and maximum expected
error parameters.

```text
    L   1.0000   0.5000   0.2500   0.1000   1.0000   0.5000   0.2500   0.1000
-----  -------  -------  -------  -------  -------  -------  -------  -------
  251     6250     3022     1305      354   21.86%   10.57%    4.56%    1.24%
  250     8894     4467     2080      637   31.10%   15.62%    7.27%    2.23%
  249     9700     4930     2291      693   33.92%   17.24%    8.01%    2.42%
  248     9854     5046     2348      712   34.46%   17.65%    8.21%    2.49%
  247     9962     5127     2390      727   34.84%   17.93%    8.36%    2.54%
[...]
    5    28595    28595    28595    28157  100.00%  100.00%  100.00%   98.47%
    4    28595    28595    28595    28461  100.00%  100.00%  100.00%   99.53%
    3    28595    28595    28595    28595  100.00%  100.00%  100.00%  100.00%
    2    28595    28595    28595    28595  100.00%  100.00%  100.00%  100.00%
    1    28595    28595    28595    28595  100.00%  100.00%  100.00%  100.00%
```


## Effect of minimum quality and length filtering

The first column indicates read lengths (Len), sorted in decreasing
order, and limited to the top *n* lengths, where n = floor(max length
/ 2) + 1. The next four columns indicate the fraction of reads that
would be retained by the `--fastq_filter` command if the reads were
truncated at length Len (option `--fastq_trunclen Len`) or at the
first position with a quality Q below 5, 10, 15 or 20 (option
`--fastq_truncqual Q`).

```text
Truncate at first Q
  Len     Q=5    Q=10    Q=15    Q=20
-----  ------  ------  ------  ------
  251   69.2%   69.2%    0.8%    0.6%
  250   93.2%   93.2%    1.6%    1.1%
  249   98.9%   98.9%    1.7%    1.2%
  248   99.8%   99.8%    1.8%    1.2%
  247   99.8%   99.8%    1.8%    1.2%
[...]
  129   99.9%   99.9%   30.3%   10.3%
  128   99.9%   99.9%   30.7%   10.4%
  127   99.9%   99.9%   31.1%   10.5%
  126   99.9%   99.9%   31.4%   10.5%
  125   99.9%   99.9%   32.1%   10.7%
```


# OPTIONS

## mandatory options

#(./fragments/option_log.md)


## core options

#(./fragments/option_fastq_ascii.md)

#(./fragments/option_fastq_qmax.md)

#(./fragments/option_fastq_qmin.md)


## secondary options

#(./fragments/option_bzip2_decompress.md)

#(./fragments/option_gzip_decompress.md)

#(./fragments/option_no_progress.md)

#(./fragments/option_quiet.md)


## ignored options

#(./fragments/option_threads_not_multithreaded.md)


# EXAMPLES

Read from *input.fastq.gz*, do not write to the
standard error (`--quiet`), and write results to *output.log*
(`--log`):

```sh
vsearch \
    --fastq_stats input.fastq.gz \
    --quiet \
    --log output.log
```

# SEE ALSO

[`vsearch-fastq_chars(1)`](./commands/vsearch-chars.1.md), 
[`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md)


#(./fragments/footer.md)
