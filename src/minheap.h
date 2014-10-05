/*
    Copyright (C) 2014 Torbjorn Rognes

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

    Contact: Torbjorn Rognes <torognes@ifi.uio.no>,
    Department of Informatics, University of Oslo,
    PO Box 1080 Blindern, NO-0316 Oslo, Norway
*/

typedef struct topscore
{
  unsigned int count;
  unsigned int seqno;
  unsigned int length;
} elem_t;


typedef struct minheap_s
{
  int alloc;
  int count;
  elem_t * array;
} minheap_t;

inline int minheap_isempty(minheap_t * m)
{
  return !m->count;
}

inline void minheap_empty(minheap_t * m)
{
  m->count = 0;
}

elem_t minheap_poplast(minheap_t * m);
void minheap_sort(minheap_t * m);
minheap_t * minheap_init(int size);
void minheap_exit(minheap_t * m);
void minheap_add(minheap_t * m, elem_t * n);
elem_t minheap_pop(minheap_t * m);

