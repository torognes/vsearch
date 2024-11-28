The *sequence* is defined as a string of IUPAC symbols
('ACGTURYSWKMDBHVN' and 'acgturyswkmdbhvn'), starting after the end of
the header line and ending before the next header line, or the file's
end. vsearch silently ignores ASCII characters 9 to 13, and exits with
an error message if ASCII characters 0 to 8, 14 to 31, '.' or '-' are
present. All other ASCII or non-ASCII characters are stripped and
complained about in a warning message.
