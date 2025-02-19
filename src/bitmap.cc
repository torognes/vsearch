/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2024, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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
#include "bitmap.h"
#include <cstring> // std::memset


auto bitmap_init(unsigned int const size) -> struct bitmap_s *
{
  auto * b = (struct bitmap_s *) xmalloc(sizeof(struct bitmap_s));
  b->size = size;
  b->bitmap = (unsigned char *) xmalloc((size + 7) / 8);
  return b;
}


auto bitmap_get(struct bitmap_s * a_bitmap, unsigned int const seed_value) -> unsigned char
{
  constexpr auto mask_111 = 7U;
  constexpr auto divider = 3U;  // divide by 8
  return (a_bitmap->bitmap[seed_value >> divider] >> (seed_value & mask_111)) & 1U;
}


auto bitmap_reset_all(struct bitmap_s * a_bitmap) -> void
{
  constexpr auto n_bits_in_a_byte = 8U;
  const auto size_in_bytes = (a_bitmap->size + n_bits_in_a_byte - 1) / n_bits_in_a_byte;
  std::memset(a_bitmap->bitmap, 0, size_in_bytes);
}


auto bitmap_set(struct bitmap_s * a_bitmap, unsigned int const seed_value) -> void
{
  constexpr auto mask_111 = 7U;
  constexpr auto divider = 3U;  // divide by 8
  a_bitmap->bitmap[seed_value >> divider] |= 1U << (seed_value & mask_111);
}


auto bitmap_free(struct bitmap_s * a_bitmap) -> void
{
  if (a_bitmap->bitmap)
    {
      xfree(a_bitmap->bitmap);
    }
  xfree(a_bitmap);
}
