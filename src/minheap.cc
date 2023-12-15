/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2021, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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
#include <cstdio>  // printf
#include <cstdlib>  // qsort()


/* implement a priority queue with a min heap binary array structure */
/* elements with the lowest count should be at the top (root) */

// refactoring: std::priority_queue (#include <queue>)

/*
  To keep track of the n best potential target sequences, we store
  them in a min heap. The root element corresponds to the least good
  target, while the best elements are found at the leaf nodes. This
  makes it simple to decide whether a new target should be included or
  not, because it just needs to be compared to the root note.  The
  list will be fully sorted before use when we want to find the best
  element and then the second best and so on.
*/

auto elem_smaller(elem_t * lhs, elem_t * rhs) -> int
{
  /* return 1 if lhs is smaller than rhs, 0 if equal or greater */
  if (lhs->count < rhs->count)
    {
      return 1;
    }
  else
    if (lhs->count > rhs->count)
      {
        return 0;
      }
    else
      if (lhs->length > rhs->length)
        {
          return 1;
        }
      else
        if (lhs->length < rhs->length)
          {
            return 0;
          }
        else
          if (lhs->seqno > rhs->seqno)
            {
              return 1;
            }
          else
            {
              return 0;
            }
}


auto minheap_compare(const void * lhs_a, const void * rhs_b) -> int
{
  auto * lhs = (elem_t *) lhs_a;
  auto * rhs = (elem_t *) rhs_b;

  /* return -1 if a is smaller than b, +1 if greater, otherwize 0 */
  /* first: lower count, larger length, lower seqno */

  if (lhs->count < rhs->count)
    {
      return -1;
    }
  else
    if (lhs->count > rhs->count)
      {
        return +1;
      }
    else
      if (lhs->length > rhs->length)
        {
          return -1;
        }
      else
        if (lhs->length < rhs->length)
          {
            return +1;
          }
        else
          if (lhs->seqno > rhs->seqno)
            {
              return -1;
            }
          else
            if (lhs->seqno < rhs->seqno)
              {
                return +1;
              }
            else
              {
                return 0;
              }
}


auto minheap_init(int size) -> minheap_t *
{
  auto * a_minheap = static_cast<minheap_t *>(xmalloc(sizeof(minheap_t)));
  a_minheap->alloc = size;
  a_minheap->array = static_cast<elem_t *>(xmalloc(size * sizeof(elem_t)));
  a_minheap->count = 0;
  return a_minheap;
}


auto minheap_exit(minheap_t * a_minheap) -> void
{
  xfree(a_minheap->array);
  xfree(a_minheap);
}


auto minheap_replaceroot(minheap_t * a_minheap, elem_t tmp) -> void
{
  /* remove the element at the root, then swap children up
     to the root and insert tmp at suitable place */

  /* start with root */
  int parent = 0;
  int nth_child = 2 * parent + 1;

  /* while at least one child */
  while (nth_child < a_minheap->count)
    {
      /* if two children: swap with the one with smallest value */
      if ((nth_child + 1 < a_minheap->count) &&
          (elem_smaller(a_minheap->array + nth_child + 1, a_minheap->array + nth_child) != 0))
        {
          ++nth_child;
        }

      /* swap parent and child if child has lower value */
      if (elem_smaller(a_minheap->array + nth_child, &tmp) != 0)
        {
          a_minheap->array[parent] = a_minheap->array[nth_child];
        }
      else
        {
          break;
        }

      /* step down */
      parent = nth_child;
      nth_child = 2 * parent + 1;
    }

  a_minheap->array[parent] = tmp;
}


auto minheap_add(minheap_t * a_minheap, elem_t * n) -> void
{
  if (a_minheap->count < a_minheap->alloc)
    {
      /* space for another item at end; swap upwards */

      int index = a_minheap->count++;
      int pos = (index - 1) / 2;
      while ((index > 0) && (elem_smaller(n, a_minheap->array + pos) != 0))
        {
          a_minheap->array[index] = a_minheap->array[pos];
          index = pos;
          pos = (index - 1) / 2;
        }
      a_minheap->array[index] = *n;
    }
  else if (elem_smaller(a_minheap->array, n) != 0)
    {
      /* replace the root if new element is larger than root */
      minheap_replaceroot(a_minheap, *n);
    }
}


auto minheap_pop(minheap_t * a_minheap) -> elem_t
{
  /* return top element and restore order */
  static const elem_t zero = {0, 0, 0};

  if (a_minheap->count != 0)
    {
      elem_t top = a_minheap->array[0];
      --a_minheap->count;
      if (a_minheap->count != 0)
        {
          const elem_t tmp = a_minheap->array[a_minheap->count];
          minheap_replaceroot(a_minheap, tmp);
        }
      return top;
    }

  return zero;
}


auto minheap_sort(minheap_t * a_minheap) -> void
{
  std::qsort(a_minheap->array, a_minheap->count, sizeof(elem_t), minheap_compare);
}


auto minheap_poplast(minheap_t * a_minheap) -> elem_t
{
  /* return top element and restore order */
  static const elem_t zero = {0, 0, 0};

  if (a_minheap->count != 0) {
    return a_minheap->array[--a_minheap->count];
  }

  return zero;
}
