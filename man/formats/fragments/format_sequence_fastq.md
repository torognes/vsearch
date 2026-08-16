The *sequence* is the line following the header line. It is a string
of IUPAC symbols ('ACGTURYSWKMDBHVN' and 'acgturyswkmdbhvn'). Carriage
returns before the newline are accepted (Windows line endings); any
other character --- including tabulations, spaces, and any non-IUPAC
letter or non-ASCII byte --- causes a fatal error. Unlike the fasta
parser, nothing is silently stripped from fastq sequences.
