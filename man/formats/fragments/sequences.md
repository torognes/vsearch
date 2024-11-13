The _sequence_ is defined as a string of IUPAC symbols
('ACGTURYSWKMDBHVN'), starting after the end of the header line and
ending before the next header line, or the file's end. vsearch
silently ignores ascii characters 9 to 13, and exits with an error
message if ascii characters 0 to 8, 14 to 31, '.' or '-' are
present. All other ascii or non-ascii characters are stripped and
complained about in a warning message.
