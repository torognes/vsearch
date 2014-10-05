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

typedef struct bitmap_s
{
  unsigned char * bitmap; /* the actual bitmap */
  unsigned int size;      /* size in bits */
} bitmap_t;

bitmap_t * bitmap_init(unsigned int size);

void bitmap_free(bitmap_t* b);

inline unsigned char bitmap_get(bitmap_t * b, unsigned int x)
{
  return (b->bitmap[x >> 3] >> (x & 7)) & 1;
}

inline void bitmap_reset_all(bitmap_t * b)
{
  memset(b->bitmap, 0, (b->size+7)/8);
}

inline void bitmap_set_all(bitmap_t * b)
{
  memset(b->bitmap, 255, (b->size+7)/8);
}

inline void bitmap_reset(bitmap_t * b, unsigned int x)
{
  b->bitmap[x >> 3] &= ~ (1 << (x & 7));
}

inline void bitmap_set(bitmap_t * b, unsigned int x)
{
  b->bitmap[x >> 3] |= 1 << (x & 7);
}

inline void bitmap_flip(bitmap_t * b, unsigned int x)
{
  b->bitmap[x >> 3] ^= 1 << (x & 7);
}
