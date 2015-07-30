/*
 * align_simd_dprofile.h
 *
 *  Created on: Jul 30, 2015
 *      Author: kaos
 */

#ifndef SRC_ALIGN_SIMD_DPROFILE_H_
#define SRC_ALIGN_SIMD_DPROFILE_H_

#include "align_simd.h"

void dprofile_fill16_aa(CELL * dprofile_word,
                        CELL * score_matrix_word,
                        BYTE * dseq);

void dprofile_fill16(CELL * dprofile_word,
                        CELL * score_matrix_word,
                        BYTE * dseq);

#endif /* SRC_ALIGN_SIMD_DPROFILE_H_ */
