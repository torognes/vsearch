/*
 Copyright (C) 2015 Jakob Frielingsdorf

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU Affero General Public License as
 published by the Free Software Foundation, either version 3 of the
 License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Affero General Public License for more details.

 You should have received a copy of the GNU Affero General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.

 Contact: Jakob Frielingsdorf <jfrielingsdorf@gmail.com>
 */

#include <stddef.h>

#include "align_simd.h"

#define MATRIX_MODE_NUC 0
#define MATRIX_MODE_AA 1

#define BLOSUM45 "blosum45"
#define BLOSUM50 "blosum50"
#define BLOSUM62 "blosum62"
#define BLOSUM80 "blosum80"
#define BLOSUM90 "blosum90"
#define PAM30 "pam30"
#define PAM70 "pam70"
#define PAM250 "pam250"

#define LINE_MAX 2048

class ScoreMatrix
{
private:
  int constant_scoring;

  void prepare_matrices(int sequence_mode);
  void finalize_matrices();

  void read_line(char line[LINE_MAX], int* symbols, char* order);

  void mat_init_from_string(const char * matrix);
  void mat_init_from_file(const char * matrix);
  void mat_init_constant_scoring(const int match, const int mismatch);
  void mat_init_buildin(const char* matrixname);

  void mat_free();

  void set16(int x, int y, long value)
  {
    score_matrix_16[(x*matrix_dim)+y] = (CELL)value;
  }

  void set64(int x, int y, long value)
  {
    score_matrix_64[(x*matrix_dim)+y] = value;
  }

public:
  int matrix_dim;

  CELL * score_matrix_16 = NULL; // short
  long * score_matrix_64 = NULL; // long

  /**
   * Assumes that amino acid sequences are aligned.
   */
  ScoreMatrix(const char* matrixname);

  /**
   * Can be used for both amino acid and nucleotide sequences.
   */
  ScoreMatrix(const int match, const int mismatch, int sequence_mode);

  ~ScoreMatrix();

  void dump_matrix();

  int is_constant_scoring()
  {
    return constant_scoring;
  }

  CELL get16(int x, int y)
  {
    return score_matrix_16[(x*matrix_dim)+y];
  }

  long get64(int x, int y)
  {
    return score_matrix_64[(x*matrix_dim)+y];
  }

};
