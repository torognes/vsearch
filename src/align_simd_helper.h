/*
 * align_simd_helper.h
 *
 *  Created on: Jul 29, 2015
 *      Author: kaos
 */

#ifndef SRC_ALIGN_SIMD_HELPER_H_
#define SRC_ALIGN_SIMD_HELPER_H_

#include <emmintrin.h>
#include <stdio.h>

#include "align_simd.h"
#include "maps.h"

#define CHANNELS 8
#define CDEPTH 4

struct s16info_s
{
  __m128i * matrix;
  __m128i * hearray;
  __m128i * dprofile;
  __m128i ** qtable;
  unsigned short * dir;
  char * qseq;
  unsigned long diralloc;

  char * cigar;
  char * cigarend;
  long cigaralloc;
  int opcount;
  char op;

  int qlen;
  int maxdlen;
  CELL penalty_gap_open_query_left;
  CELL penalty_gap_open_target_left;
  CELL penalty_gap_open_query_interior;
  CELL penalty_gap_open_target_interior;
  CELL penalty_gap_open_query_right;
  CELL penalty_gap_open_target_right;
  CELL penalty_gap_extension_query_left;
  CELL penalty_gap_extension_target_left;
  CELL penalty_gap_extension_query_interior;
  CELL penalty_gap_extension_target_interior;
  CELL penalty_gap_extension_query_right;
  CELL penalty_gap_extension_target_right;
};

//void _mm_print(__m128i x)
//{
//  unsigned short * y = (unsigned short*)&x;
//  for (int i = 0; i<8; i++)
//    printf("%s%6d", (i>0 ? " " : ""), y[7-i]);
//}
//
//void _mm_print16(__m128i x)
//{
//  unsigned char * y = (unsigned char*)&x;
//  for (int i = 0; i<16; i++)
//    printf("%s%2d", (i>0 ? " " : ""), y[15-i]);
//}
//
//void _mm_print2(__m128i x)
//{
//  signed short * y = (signed short*)&x;
//  for (int i = 0; i<8; i++)
//    printf("%s%2d", (i>0 ? " " : ""), y[7-i]);
//}
//
//void dprofile_dump16(CELL * dprofile)
//{
//  char * s = sym_nt_4bit;
//  printf("\ndprofile:\n");
//  for (int i = 0; i<16; i++)
//    {
//      printf("%c: ", s[i]);
//      for (int k = 0; k<CDEPTH; k++)
//        {
//          printf("[");
//          for (int j = 0; j<CHANNELS; j++)
//            printf(" %3d", dprofile[CHANNELS*CDEPTH*i+CHANNELS*k+j]);
//          printf("]");
//        }
//      printf("\n");
//    }
//}

#endif /* SRC_ALIGN_SIMD_HELPER_H_ */
