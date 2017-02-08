/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2017, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
  All rights reserved.

  Contact: Torbjorn Rognes <torognes@ifi.uio.no>,
  Department of Informatics, University of Oslo,
  PO Box 1080 Blindern, NO-0316 Oslo, Norway

  This software is dual-licensed and available under a choice
  of one of two licenses, either under the terms of the GNU
  General Public License version 3 or the BSD 2-Clause License.


  GNU General Public License version 3

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.


  The BSD 2-Clause License

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions
  are met:

  1. Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

  2. Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
  COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
  POSSIBILITY OF SUCH DAMAGE.

*/

#include "vsearch.h"

/* implement a priority queue with a min heap binary array structure */
/* elements with the lowest count should be at the top (root) */

/*
  To keep track of the n best potential target sequences, we store
  them in a min heap. The root element corresponds to the least good
  target, while the best elements are found at the leaf nodes. This
  makes it simple to decide whether a new target should be included or
  not, because it just needs to be compared to the root note.  The
  list will be fully sorted before use when we want to find the best
  element and then the second best and so on.
*/

int
elem_smaller(elem_t * a, elem_t * b)
{
  /* return 1 if a is smaller than b, 0 if equal or greater */
  if (a->count < b->count)
    return 1;
  else
    if (a->count > b->count)
      return 0;
    else
      if (a->length > b->length)
        return 1;
      else
        if (a->length < b->length)
          return 0;
        else
          if (a->seqno > b->seqno)
            return 1;
          else
            return 0;
}

int minheap_compare(const void * a, const void * b)
{
  elem_t * x = (elem_t*) a;
  elem_t * y = (elem_t*) b;

  /* return -1 if a is smaller than b, +1 if greater, otherwize 0 */
  /* first: lower count, larger length, lower seqno */

  if (x->count < y->count)
    return -1;
  else
    if (x->count > y->count)
      return +1;
    else
      if (x->length > y->length)
        return -1;
      else
        if (x->length < y->length)
          return +1;
        else
          if (x->seqno > y->seqno)
            return -1;
          else
            if (x->seqno < y->seqno)
              return +1;
            else
              return 0;
}

minheap_t *
minheap_init(int size)
{
  minheap_t * m = (minheap_t *) xmalloc(sizeof(minheap_t));
  m->alloc = size;
  m->array = (elem_t *) xmalloc(size * sizeof(elem_t));
  m->count = 0;
  return m;
}

void
minheap_exit(minheap_t * m)
{
  xfree(m->array);
  xfree(m);
}

static int swaps = 0;

void
minheap_replaceroot(minheap_t * m, elem_t tmp)
{
  /* remove the element at the root, then swap children up
     to the root and insert tmp at suitable place */

  /* start with root */
  int p = 0;
  int c = 2*p+1;
  
  /* while at least one child */
  while (c < m->count)
    {
      /* if two children: swap with the one with smallest value */
      if ((c + 1 < m->count) && (elem_smaller(m->array + c + 1, m->array + c)))
        c++;
      
      /* swap parent and child if child has lower value */
      if (elem_smaller(m->array + c, &tmp))
        {
          m->array[p] = m->array[c];
          swaps++;
        }
      else
        break;
      
      /* step down */
      p = c;
      c = 2*p+1;
    }
  
  m->array[p] = tmp;
}


void
minheap_add(minheap_t * m, elem_t * n)
{
  if (m->count < m->alloc)
    {
      /* space for another item at end; swap upwards */

      int i = m->count++;
      int p = (i-1)/2;
      while ((i>0) && elem_smaller(n, m->array+p))
        {
          m->array[i] = m->array[p];
          i = p;
          p = (i-1)/2;
          swaps++;
        }
      m->array[i] = *n;
    }
  else if (elem_smaller(m->array, n))
    {
      /* replace the root if new element is larger than root */
      minheap_replaceroot(m, *n);
    }
}

#if 0
inline int
minheap_isempty(minheap_t * m)
{
  return !m->count;
}

inline void
minheap_empty(minheap_t * m)
{
  m->count = 0;
}
#endif

elem_t
minheap_pop(minheap_t * m)
{
  /* return top element and restore order */
  static elem_t zero = {0,0,0};

  if (m->count)
    {
      elem_t top = m->array[0];
      m->count--;
      if (m->count)
        {
          elem_t tmp = m->array[m->count];
          minheap_replaceroot(m, tmp);
        }
      return top;
    }
  else
    return zero;
}

void
minheap_sort(minheap_t * m)
{
  qsort(m->array, m->count, sizeof(elem_t), minheap_compare);
}

void
minheap_dump(minheap_t * m)
{
  for(int i=0; i < m->count; i++)
    printf("%s%u", i>0 ? " " : "", m->array[i].count);
  printf("\n");
}

elem_t
minheap_poplast(minheap_t * m)
{
  /* return top element and restore order */
  static elem_t zero = {0,0,0};

  if (m->count)
    return m->array[--m->count];
  else
    return zero;
}


void
minheap_test()
{
  minheap_t * m = minheap_init(10000000);

  int samples = 10000000;
  swaps = 0;
  
  for(int i=samples; i>=0; i--)
    {
      elem_t x = {(unsigned int)(rand()),0,1};
      minheap_add(m, & x);
    }

  minheap_sort(m);

  while(! minheap_isempty(m))
    {
      elem_t x = minheap_poplast(m);
      printf("%u\n", x.count);
    }

  printf("Swaps: %d\n\n", swaps);

  minheap_exit(m);
}
