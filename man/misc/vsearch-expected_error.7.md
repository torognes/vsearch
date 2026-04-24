% vsearch-expected_error(7) version 2.31.0 | vsearch file formats
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(../commands/fragments/date.md)

# NAME

expected error --- a quality summary metric for fastq sequences


# DESCRIPTION

The *expected error* (EE) of a fastq read is a single number that
summarises the overall quality of that read. It is defined as the sum
of the per-base error probabilities across the entire sequence:

```text
EE = P_e(1) + P_e(2) + ... + P_e(L)
```

where *L* is the length of the read and *P_e(i)* is the probability
that base *i* was incorrectly sequenced.


## Computing error probabilities from quality scores

Each base quality score *Q* stored in a fastq file encodes its error
probability according to the Phred scale:

```text
P_e = 10^(-Q / 10)
```

For example:

| Q score | Error probability | Accuracy |
|---------|-------------------|----------|
| Q10     | 0.1               | 90%      |
| Q20     | 0.01              | 99%      |
| Q30     | 0.001             | 99.9%    |
| Q40     | 0.0001            | 99.99%   |

Since error probabilities are always positive, the expected error of a
sequence is always greater than zero. It is at most equal to the
sequence length, which occurs when every base has an error probability
of 1.0.


## Why expected error is preferred over average quality

Averaging Phred quality scores directly is mathematically incorrect:
quality scores are on a logarithmic scale, while error probabilities are
linear. The expected error sums all per-base error probabilities, giving
each base equal weight and producing a meaningful global quality
estimate.


## Poisson interpretation

Because sequencing errors are approximately independent, the expected
error can be used as the *lambda* (λ) parameter of the Poisson
distribution to estimate the probability of observing exactly *k*
errors in a read. For a read with EE = 1.0:

- 36.8% chance of zero errors,
- 36.8% chance of one error,
- 18.4% chance of two errors,
- 6.1% chance of three errors,
- 1.5% chance of four errors,
- 0.3% chance of five errors,
- etc.


## vsearch options using expected error

`--fastq_maxee` *real*
: Discard reads whose total expected error exceeds *real*. Applied
  after any trimming step.

`--fastq_maxee_rate` *real*
: Discard reads whose average expected error per base (EE divided by
  read length) exceeds *real* (values from 0.0 to 1.0).

`--fastq_truncee` *real*
: Truncate reads at the first position where the cumulative expected
  error would exceed *real*.

`--fastq_truncee_rate` *real*
: Truncate reads at the first position where the average expected
  error per base (cumulative EE divided by current length) would
  exceed *real*.

`--fastq_eeout` / `--eeout`
: Annotate output sequence headers with the expected error value as
  `;ee=float`. Use `--xee` to remove this annotation from headers.


# SEE ALSO

[`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md),
[`vsearch-fastq_eestats(1)`](../commands/vsearch-fastq_eestats.1.md),
[`vsearch-fastq_eestats2(1)`](../commands/vsearch-fastq_eestats2.1.md),
[`vsearch-fastq_filter(1)`](../commands/vsearch-fastq_filter.1.md),
[`vsearch-fastq_mergepairs(1)`](../commands/vsearch-fastq_mergepairs.1.md),
[`vsearch-fastx_filter(1)`](../commands/vsearch-fastx_filter.1.md)


#(../commands/fragments/footer.md)
