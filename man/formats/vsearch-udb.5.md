% vsearch-udb(5) version 2.30.4 | vsearch file formats
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(../commands/fragments/date.md)

# NAME

udb --- a binary file format containing fasta sequences and a k-mer index


# DESCRIPTION

The UDB (USEARCH database) format is a binary file format used to store
fasta sequences together with a pre-computed k-mer index. Building a UDB
file from a fasta reference database with `--makeudb_usearch` allows
subsequent search commands (such as `--usearch_global`) to load the index
directly, avoiding the cost of recomputing it at each run.

UDB files are not portable across machines with different native byte
orders. All multi-byte numeric values are stored in *little-endian* byte
order (the least significant byte first), which is the native byte order
of x86 and ARM processors. This contrasts with *big-endian* order (most
significant byte first), which is used for example by the SFF format (see
[`vsearch-sff(5)`](./vsearch-sff.5.md)).

UDB files cannot be read from pipes; a seekable file path must be
provided.

The file consists of nine sequential sections described below.


## Section 1 --- Main Header

The main header is exactly 200 bytes long (50 × `uint32_t`):

```text
 buffer[0]   uint32_t   magic number: 0x55444246 ("UDBF")
 buffer[1]   uint32_t   reserved, set to 0
 buffer[2]   uint32_t   sequence index bits: 32
 buffer[3]   uint32_t   reserved, set to 0
 buffer[4]   uint32_t   word length (k-mer size, default 8)
 buffer[5]   uint32_t   dbstep: 1
 buffer[6]   uint32_t   dbaccelpct: 100
 buffer[7]   uint32_t   reserved, set to 0
 ...
 buffer[11]  uint32_t   slots: 0
 buffer[12]  uint32_t   reserved, set to 0
 buffer[13]  uint32_t   number of sequences
 buffer[14]  uint32_t   reserved, set to 0
 ...
 buffer[17]  uint32_t   alphabet: 0x0000746e ("nt")
 buffer[18]  uint32_t   reserved, set to 0
 ...
 buffer[49]  uint32_t   end marker: 0x55444266 ("UDBf")
```

All fields not listed explicitly are reserved and set to zero.

- The `magic number` field value is 0x55444246, the little-endian
  encoding of the ASCII string "UDBF".
- The `word length` field stores the k-mer size used to build the index
  (between 3 and 15, default 8; controllable with `--wordlength`).
- The `dbstep` field is always 1 and the `dbaccelpct` field is always
  100, reflecting that every sequence is indexed with full coverage.
- The `number of sequences` field gives the total number of sequences
  stored in the file.
- The `alphabet` field value is 0x0000746e, the little-endian encoding
  of the ASCII string "nt" (nucleotide), indicating that this is a
  nucleotide database.
- The `end marker` field value is 0x55444266, the little-endian encoding
  of the ASCII string "UDBf", marking the end of the main header.


## Section 2 --- Word Match Counts

This section contains `4^wordlength` consecutive `uint32_t` values
(one per possible k-mer):

```text
 kmercount[0]                uint32_t
 kmercount[1]                uint32_t
 ...
 kmercount[4^wordlength - 1] uint32_t
```

Each value `kmercount[i]` stores the number of sequences in the database
that contain k-mer `i` at least once. With the default word length of 8,
this section contains 65,536 values (256 KB).


## Section 3 --- UDB3 Marker

A single `uint32_t` sentinel value:

```text
 marker   uint32_t   0x55444233 ("UDB3")
```

The value 0x55444233 is the little-endian encoding of the ASCII string
"UDB3".


## Section 4 --- Word Index

This section stores, for each k-mer, the sorted list of 0-based sequence
numbers of all sequences containing that k-mer. The lists are
concatenated in k-mer order, with no delimiters:

```text
 seqno[0]   uint32_t   first sequence number for k-mer 0
 seqno[1]   uint32_t   second sequence number for k-mer 0
 ...        (kmercount[0] entries for k-mer 0)
 seqno[.]   uint32_t   first sequence number for k-mer 1
 ...        (kmercount[1] entries for k-mer 1)
 ...
```

The total number of `uint32_t` values in this section equals the sum of
all `kmercount` values from section 2. A k-mer with a count of zero
contributes no bytes to this section. Sequence numbers use 0-based
indexing.


## Section 5 --- UDB4 Header

A fixed-size block of 32 bytes (8 × `uint32_t`):

```text
 buffer[0]   uint32_t   magic number: 0x55444234 ("UDB4")
 buffer[1]   uint32_t   constant: 0x005e0db3
 buffer[2]   uint32_t   number of sequences
 buffer[3]   uint32_t   total nucleotide count (low 32 bits)
 buffer[4]   uint32_t   total nucleotide count (high 32 bits)
 buffer[5]   uint32_t   total header characters (low 32 bits)
 buffer[6]   uint32_t   total header characters (high 32 bits)
 buffer[7]   uint32_t   constant: 0x005e0db4
```

- The `magic number` field value is 0x55444234, the little-endian
  encoding of the ASCII string "UDB4".
- The `number of sequences` field duplicates the value from the main
  header.
- The `total nucleotide count` is a `uint64_t` stored as two consecutive
  `uint32_t` values (low word first). It gives the total number of bases
  across all sequences, which equals the size of section 9.
- The `total header characters` is a `uint64_t` stored as two consecutive
  `uint32_t` values (low word first). It gives the total number of bytes
  in section 7, including one null terminator per sequence.


## Section 6 --- Header Index

This section contains `seqcount` consecutive `uint32_t` values:

```text
 offset[0]            uint32_t   byte offset of sequence 0 header
 offset[1]            uint32_t   byte offset of sequence 1 header
 ...
 offset[seqcount-1]   uint32_t   byte offset of last sequence header
```

Each value is the byte offset of the corresponding sequence's header
string within section 7. Offsets use 0-based indexing relative to the
start of section 7.


## Section 7 --- Headers

This section contains `header_characters` bytes: the ASCII sequence
headers concatenated in order, each terminated by a null byte (`\0`).
The section is not padded. There are no '>' characters; the header
strings correspond to the part of fasta header lines after the leading
'>'.


## Section 8 --- Sequence Lengths

This section contains `seqcount` consecutive `uint32_t` values:

```text
 length[0]            uint32_t   length of sequence 0
 length[1]            uint32_t   length of sequence 1
 ...
 length[seqcount-1]   uint32_t   length of last sequence
```

Each value is the number of bases in the corresponding sequence.


## Section 9 --- Sequences

This section contains `ntcount` bytes: the ASCII nucleotide sequences
concatenated in order. The section is not null-terminated and not padded.
Sequences can contain uppercase and lowercase letters (soft masking is
preserved). T and U are treated as equivalent by vsearch.


# EXAMPLES

Build a UDB file from a fasta reference database:

```sh
vsearch \
    --makeudb_usearch db.fasta \
    --output db.udb
```

Use a UDB file to search query sequences against the indexed database:

```sh
vsearch \
    --usearch_global queries.fasta \
    --db db.udb \
    --id 0.97 \
    --blast6out results.tsv
```

Inspect the content of a UDB file:

```sh
vsearch --udbinfo db.udb
```

Extract the sequences stored in a UDB file back to fasta:

```sh
vsearch --udb2fasta db.udb --output db_extracted.fasta
```


# SEE ALSO

[`vsearch-makeudb_usearch(1)`](../commands/vsearch-makeudb_usearch.1.md),
[`vsearch-udb2fasta(1)`](../commands/vsearch-udb2fasta.1.md),
[`vsearch-udbinfo(1)`](../commands/vsearch-udbinfo.1.md),
[`vsearch-udbstats(1)`](../commands/vsearch-udbstats.1.md),
[`vsearch-fasta(5)`](./vsearch-fasta.5.md),
[`vsearch-sff(5)`](./vsearch-sff.5.md)


#(../commands/fragments/footer.md)
