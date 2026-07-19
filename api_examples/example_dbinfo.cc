/*
 * example_dbinfo.cc — Database query and indexing library surface.
 *
 * The global sequence database exposes a read API (counts, per-record
 * accessors, quality) plus preparation primitives (sort variants, incremental
 * k-mer indexing) that a library embedder drives directly but the CLI only
 * uses internally. None of these accessors are checked by the other examples,
 * which read headers/sequences but never verify the length/nucleotide/
 * abundance/quality bookkeeping. Every check here is self-validating against
 * values known from the input, not against native vsearch output.
 *
 * Two of these checks are regression guards for library-only bugs fixed
 * alongside them (db.cc):
 *   - db_read() with a default Parameters used to discard every sequence: the
 *     -1 "unset" opt_minseqlength sentinel was cast to size_t (SIZE_MAX), so
 *     every length compared as "too short". db_read now treats a non-positive
 *     minimum as "no lower bound" (test_db_read_fasta_accessors uses a default
 *     Parameters and expects all sequences to load).
 *   - db_add(is_fastq=true, ...) stored the quality string but left the global
 *     is_fastq flag false, so db_is_fastq()/db_getquality() could not reach it.
 *     db_add now sets the flag (test_db_add_fastq_quality).
 *
 * Build:  g++ -std=c++11 -O3 -I../src -o example_dbinfo example_dbinfo.cc ../src/libvsearch.a -lpthread -ldl
 * Run:    ./example_dbinfo
 */

#include "vsearch_api.h"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>


/* --- Test 1: db_read (FASTA) + statistical accessors ---
   chimera_ref.fasta holds 6 sequences of 300 nt each with headers "ref1".."ref6"
   (4 characters). Verify every summary accessor against those known values.
   Also a regression guard: this uses a *default* Parameters (opt_minseqlength
   left at its -1 sentinel), so all 6 sequences must load — previously db_read
   discarded every sequence in this case. */
static int test_db_read_fasta_accessors()
{
  int failures = 0;

  struct Parameters parameters;
  VsearchSession const session(parameters);

  Database db;
  db.read("data/chimera_ref.fasta", 0, parameters);

  struct { char const * name; uint64_t got; uint64_t want; } const checks[] = {
    { "sequencecount",  db.getsequencecount(),  6    },
    { "nucleotidecount", db.getnucleotidecount(), 1800 },
    { "longestsequence", db.getlongestsequence(), 300  },
    { "shortestsequence", db.getshortestsequence(), 300 },
    { "longestheader",   db.getlongestheader(),   4    },
  };
  for (auto const & check : checks)
    {
      if (check.got != check.want)
        {
          std::fprintf(stderr, "FAIL: db_get%s() = %lu, expected %lu\n",
                       check.name, (unsigned long) check.got, (unsigned long) check.want);
          ++failures;
        }
    }

  /* Per-record accessors: header length, sequence length, default abundance. */
  for (uint64_t seqno = 0; seqno < db.getsequencecount(); ++seqno)
    {
      if (db.getsequencelen(seqno) != 300)
        {
          std::fprintf(stderr, "FAIL: db_getsequencelen(%lu) = %lu, expected 300\n",
                       (unsigned long) seqno, (unsigned long) db.getsequencelen(seqno));
          ++failures;
        }
      if (db.getheaderlen(seqno) != std::strlen(db.getheader(seqno)))
        {
          std::fprintf(stderr, "FAIL: db_getheaderlen(%lu) disagrees with strlen(header)\n",
                       (unsigned long) seqno);
          ++failures;
        }
      if (db.getabundance(seqno) != 1)
        {
          std::fprintf(stderr, "FAIL: db_getabundance(%lu) = %lu, expected 1 (no ;size=)\n",
                       (unsigned long) seqno, (unsigned long) db.getabundance(seqno));
          ++failures;
        }
    }

  if (db.is_fastq())
    {
      std::fprintf(stderr, "FAIL: db_is_fastq() true after reading a FASTA file\n");
      ++failures;
    }

  db.clear();

  if (failures == 0)
    {
      std::fprintf(stderr, "PASS: db_read(FASTA) + accessors report the known statistics\n");
    }
  return failures;
}


/* --- Test 2: db_read (FASTQ) exposes quality via db_getquality ---
   A FASTQ file sets db_is_fastq() and makes db_getquality() return the stored
   quality string, whose length matches the sequence. */
static int test_db_read_fastq_quality()
{
  int failures = 0;

  struct Parameters parameters;
  VsearchSession const session(parameters);

  Database db;
  db.read("data/merge_fwd.fastq", 0, parameters);

  if (not db.is_fastq())
    {
      std::fprintf(stderr, "FAIL: db_is_fastq() false after reading a FASTQ file\n");
      ++failures;
    }
  if (db.getsequencecount() < 1)
    {
      std::fprintf(stderr, "FAIL: FASTQ read produced no sequences\n");
      ++failures;
    }
  else
    {
      char const * const quality = db.getquality(0);
      if (quality == nullptr)
        {
          std::fprintf(stderr, "FAIL: db_getquality(0) is null for a FASTQ database\n");
          ++failures;
        }
      else if (std::strlen(quality) != db.getsequencelen(0))
        {
          std::fprintf(stderr, "FAIL: quality length %lu != sequence length %lu\n",
                       (unsigned long) std::strlen(quality),
                       (unsigned long) db.getsequencelen(0));
          ++failures;
        }
    }

  db.clear();

  if (failures == 0)
    {
      std::fprintf(stderr, "PASS: db_read(FASTQ) sets db_is_fastq() and exposes quality\n");
    }
  return failures;
}


/* --- Test 2b: db_add(is_fastq=true) exposes quality (regression guard) ---
   Building a FASTQ database directly with db_add() must set db_is_fastq() and
   make db_getquality() return the stored quality. Previously db_add stored the
   quality but left the global flag false, so the quality was unreachable. */
static int test_db_add_fastq_quality()
{
  int failures = 0;

  struct Parameters parameters;
  VsearchSession const session(parameters);

  Database db;
  db.init();
  char const * const header = "read1";
  char const * const sequence = "ACGTACGTACGT";
  char const * const quality  = "IIIIIIIIIIII";
  db.add(true, SeqRecord{View<char>{header, std::strlen(header)},
                         View<char>{sequence, std::strlen(sequence)},
                         View<char>{quality, std::strlen(sequence)}}, 1);

  if (not db.is_fastq())
    {
      std::fprintf(stderr, "FAIL: db_is_fastq() false after db_add(is_fastq=true)\n");
      ++failures;
    }
  char const * const stored = db.getquality(0);
  if (stored == nullptr)
    {
      std::fprintf(stderr, "FAIL: db_getquality(0) null after db_add(is_fastq=true)\n");
      ++failures;
    }
  else if (std::strcmp(stored, quality) != 0)
    {
      std::fprintf(stderr, "FAIL: db_getquality(0) = '%s', expected '%s'\n", stored, quality);
      ++failures;
    }

  db.clear();

  if (failures == 0)
    {
      std::fprintf(stderr, "PASS: db_add(is_fastq=true) sets db_is_fastq() and exposes quality\n");
    }
  return failures;
}


/* Six synthetic records with distinct lengths and abundances. */
struct record_s { std::string header; std::string sequence; int64_t abundance; };

static std::vector<record_s> make_records()
{
  int const lengths[]     = { 40, 20, 60, 30, 50, 25 };
  int64_t const abunds[]  = {  5, 50,  3, 20, 100, 8 };
  std::vector<record_s> records;
  for (int idx = 0; idx < 6; ++idx)
    {
      record_s record;
      record.header = "seq" + std::to_string(idx);
      /* deterministic non-repetitive-ish sequence of the requested length */
      static char const bases[] = "ACGT";
      for (int pos = 0; pos < lengths[idx]; ++pos)
        {
          record.sequence += bases[(pos + (3 * idx)) % 4];
        }
      record.abundance = abunds[idx];
      records.push_back(record);
    }
  return records;
}

static void load_records(Database & db, std::vector<record_s> const & records)
{
  db.init();
  for (auto const & record : records)
    {
      db.add(false, SeqRecord{View<char>{record.header.c_str(), record.header.size()},
                              View<char>{record.sequence.c_str(), record.sequence.size()},
                              View<char>{}}, record.abundance);
    }
}


/* --- Test 3: db_add bookkeeping (abundance, per-record lengths, totals) --- */
static int test_db_add_accessors()
{
  int failures = 0;

  struct Parameters parameters;
  VsearchSession const session(parameters);

  std::vector<record_s> const records = make_records();
  Database db;
  load_records(db, records);

  uint64_t expected_nt = 0;
  uint64_t expected_longest = 0;
  uint64_t expected_shortest = UINT64_MAX;
  for (auto const & record : records)
    {
      expected_nt += record.sequence.size();
      expected_longest = std::max<uint64_t>(expected_longest, record.sequence.size());
      expected_shortest = std::min<uint64_t>(expected_shortest, record.sequence.size());
    }

  if (db.getsequencecount() != records.size())
    {
      std::fprintf(stderr, "FAIL: db_getsequencecount() = %lu, expected %zu\n",
                   (unsigned long) db.getsequencecount(), records.size());
      ++failures;
    }
  if (db.getnucleotidecount() != expected_nt)
    {
      std::fprintf(stderr, "FAIL: db_getnucleotidecount() = %lu, expected %lu\n",
                   (unsigned long) db.getnucleotidecount(), (unsigned long) expected_nt);
      ++failures;
    }
  if (db.getlongestsequence() != expected_longest ||
      db.getshortestsequence() != expected_shortest)
    {
      std::fprintf(stderr, "FAIL: longest/shortest = %lu/%lu, expected %lu/%lu\n",
                   (unsigned long) db.getlongestsequence(),
                   (unsigned long) db.getshortestsequence(),
                   (unsigned long) expected_longest, (unsigned long) expected_shortest);
      ++failures;
    }
  for (uint64_t seqno = 0; seqno < db.getsequencecount(); ++seqno)
    {
      if (db.getabundance(seqno) != static_cast<uint64_t>(records[seqno].abundance) ||
          db.getsequencelen(seqno) != records[seqno].sequence.size())
        {
          std::fprintf(stderr, "FAIL: record %lu abundance/len mismatch\n",
                       (unsigned long) seqno);
          ++failures;
        }
    }

  db.clear();

  if (failures == 0)
    {
      std::fprintf(stderr, "PASS: db_add() tracks abundance, lengths and totals\n");
    }
  return failures;
}


/* Collect the multiset of sequence lengths currently in the database. */
static std::vector<uint64_t> current_lengths(Database const & db)
{
  std::vector<uint64_t> lengths;
  for (uint64_t seqno = 0; seqno < db.getsequencecount(); ++seqno)
    {
      lengths.push_back(db.getsequencelen(seqno));
    }
  std::sort(lengths.begin(), lengths.end());
  return lengths;
}


/* --- Test 4: sort ordering contracts (and multiset preservation) ---
   db_sortbyabundance / db_sortbylength / db_sortbylength_shortest_first each
   reorder the database. Verify the documented ordering after each and that no
   sequence is lost or duplicated (the multiset of lengths is invariant). */
static int test_sort_contracts()
{
  int failures = 0;

  struct Parameters parameters;
  VsearchSession const session(parameters);

  std::vector<record_s> const records = make_records();
  Database db;
  load_records(db, records);
  std::vector<uint64_t> const reference_lengths = current_lengths(db);

  /* db_sortbyabundance: abundance non-increasing */
  db.sortbyabundance(parameters);
  for (uint64_t seqno = 1; seqno < db.getsequencecount(); ++seqno)
    {
      if (db.getabundance(seqno) > db.getabundance(seqno - 1))
        {
          std::fprintf(stderr, "FAIL: db_sortbyabundance not non-increasing at %lu\n",
                       (unsigned long) seqno);
          ++failures;
        }
    }

  /* db_sortbylength: length non-increasing */
  db.sortbylength(parameters);
  for (uint64_t seqno = 1; seqno < db.getsequencecount(); ++seqno)
    {
      if (db.getsequencelen(seqno) > db.getsequencelen(seqno - 1))
        {
          std::fprintf(stderr, "FAIL: db_sortbylength not non-increasing at %lu\n",
                       (unsigned long) seqno);
          ++failures;
        }
    }

  /* db_sortbylength_shortest_first: length non-decreasing */
  db.sortbylength_shortest_first(parameters);
  for (uint64_t seqno = 1; seqno < db.getsequencecount(); ++seqno)
    {
      if (db.getsequencelen(seqno) < db.getsequencelen(seqno - 1))
        {
          std::fprintf(stderr, "FAIL: db_sortbylength_shortest_first not non-decreasing at %lu\n",
                       (unsigned long) seqno);
          ++failures;
        }
    }

  if (current_lengths(db) != reference_lengths)
    {
      std::fprintf(stderr, "FAIL: sorting changed the multiset of sequences\n");
      ++failures;
    }

  db.clear();

  if (failures == 0)
    {
      std::fprintf(stderr, "PASS: sort variants honour their ordering and preserve the sequences\n");
    }
  return failures;
}


/* Run one query against the currently-indexed database and return the sorted
   list of "target:%.2f id" hit descriptors, so two index builds can be compared. */
static std::vector<std::string> search_hits(struct Parameters const & parameters,
                                            struct Dbindex const & dbindex,
                                            struct Database const & db,
                                            char const * query)
{
  struct search_session_s * const session = search_session_alloc();
  search_session_init(session, parameters, dbindex, db);

  struct search_result_s results[16];
  int count = 0;
  search_session_single(session, query, "q",
                        static_cast<int>(std::strlen(query)), 1,
                        results, 16, &count);

  std::vector<std::string> hits;
  for (int idx = 0; idx < count; ++idx)
    {
      char buffer[64];
      std::snprintf(buffer, sizeof(buffer), "%d:%.2f", results[idx].target, results[idx].id);
      hits.push_back(buffer);
    }
  std::sort(hits.begin(), hits.end());

  search_session_cleanup(session);
  search_session_free(session);
  return hits;
}


/* --- Test 5: dbindex_addsequence (incremental) == dbindex_addallsequences ---
   dbindex_addsequence is the incremental indexing primitive (used internally by
   clustering and de novo chimera detection) that no other example exercises.
   Build the index both ways over the same database and require identical search
   results. */
static int test_incremental_indexing()
{
  int failures = 0;

  struct Parameters parameters;
  parameters.opt_wordlength = 8;
  parameters.opt_id = 0.5;
  parameters.opt_maxaccepts = 8;
  parameters.opt_maxrejects = 32;
  VsearchSession const session(parameters);

  std::vector<std::string> ref_labels, ref_seqs;
  {
    std::FILE * const input = std::fopen("data/chimera_ref.fasta", "r");
    if (input == nullptr)
      {
        std::fprintf(stderr, "FAIL: cannot open data/chimera_ref.fasta\n");
        return 1;
      }
    char line[65536];
    std::string label, seq;
    while (std::fgets(line, sizeof(line), input) != nullptr)
      {
        char * newline = std::strchr(line, '\n');
        if (newline != nullptr) { *newline = '\0'; }
        if (line[0] == '>')
          {
            if (not label.empty()) { ref_labels.push_back(label); ref_seqs.push_back(seq); }
            label = line + 1; seq.clear();
          }
        else { seq += line; }
      }
    if (not label.empty()) { ref_labels.push_back(label); ref_seqs.push_back(seq); }
    std::fclose(input);
  }

  Database db;
  auto load_db = [&]() {
    db.init();
    for (size_t idx = 0; idx < ref_labels.size(); ++idx)
      {
        db.add(false, SeqRecord{View<char>{ref_labels[idx].c_str(), ref_labels[idx].size()},
                                View<char>{ref_seqs[idx].c_str(), ref_seqs[idx].size()},
                                View<char>{}}, 1);
      }
    dust_all(db, parameters);
  };

  char const * const query = ref_seqs[0].c_str();

  /* (a) batch indexing */
  load_db();
  Dbindex dbindex;
  dbindex.prepare(1, parameters.opt_dbmask, db, parameters);
  dbindex.add_all_sequences(parameters.opt_dbmask, db, parameters);
  std::vector<std::string> const hits_batch = search_hits(parameters, dbindex, db, query);
  dbindex.clear();

  /* (b) incremental indexing: one dbindex_addsequence() per sequence */
  dbindex.prepare(1, parameters.opt_dbmask, db, parameters);
  for (uint64_t seqno = 0; seqno < db.getsequencecount(); ++seqno)
    {
      dbindex.add_sequence(static_cast<unsigned int>(seqno), parameters.opt_dbmask, db);
    }
  std::vector<std::string> const hits_incremental = search_hits(parameters, dbindex, db, query);
  dbindex.clear();

  db.clear();

  if (hits_batch.empty())
    {
      std::fprintf(stderr, "FAIL: batch index produced no hits (test is vacuous)\n");
      ++failures;
    }
  if (hits_batch != hits_incremental)
    {
      std::fprintf(stderr, "FAIL: incremental indexing gave different hits than batch\n");
      std::fprintf(stderr, "      batch=%zu incremental=%zu\n",
                   hits_batch.size(), hits_incremental.size());
      ++failures;
    }

  if (failures == 0)
    {
      std::fprintf(stderr, "PASS: incremental dbindex_addsequence matches dbindex_addallsequences (%zu hits)\n",
                   hits_batch.size());
    }
  return failures;
}


int main()
{
  int failures = 0;
  failures += test_db_read_fasta_accessors();
  failures += test_db_read_fastq_quality();
  failures += test_db_add_fastq_quality();
  failures += test_db_add_accessors();
  failures += test_sort_contracts();
  failures += test_incremental_indexing();
  return failures == 0 ? 0 : 1;
}
