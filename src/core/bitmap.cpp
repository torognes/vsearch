/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2026, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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

#include "core/bitmap.hpp"
#include "utils/span.hpp"  // Span
#include <algorithm>  // std::fill
#include <cassert>
#include <cstddef>  // std::ptrdiff_t
#include <cstdint>  // uint64_t
#include <cstring>  // std::memcpy
#include <iterator>  // std::next


// anonymous namespace: limit visibility and usage to this translation unit
namespace {

  /* Number of trailing zero bits, i.e. the index of the lowest set bit.
     Undefined for zero, as the standard function it stands in for is; the one
     caller tests the word first.

     C++20 refactoring: replace with std::countr_zero (<bit>).

     The guarded-builtin shape is the one utils/os_byteswap.cpp already uses for
     bswap_*: the builtin where the compiler provides it -- every compiler
     vsearch is built with, mingw included -- and a portable loop otherwise, so
     the file does not depend on a GNU extension being present. */
  auto countr_zero(uint64_t const word) noexcept -> unsigned int
  {
#if defined(__GNUC__) || defined(__clang__)
    return static_cast<unsigned int>(__builtin_ctzll(word));
#else
    auto count = 0U;
    while (((word >> count) & 1ULL) == 0ULL)
      {
        ++count;
      }
    return count;
#endif
  }

}  // end of anonymous namespace


Bitmap::Bitmap(unsigned int const size)
{
  constexpr auto divider = 8U;
  constexpr auto padding = divider - 1U;
  bitmap_.assign((size + padding) / divider, 0U);
}


auto Bitmap::empty() const -> bool
{
  return bitmap_.empty();
}


auto Bitmap::data() const -> unsigned char const *
{
  return bitmap_.data();
}


auto Bitmap::is_set(unsigned int const seed_value) const -> bool
{
  constexpr auto mask_111 = 7U;
  constexpr auto divider = 3U;  // divide by 8
  return ((bitmap_[seed_value >> divider] >> (seed_value & mask_111)) & 1U) != 0U;
}


auto Bitmap::reset_all() -> void
{
  std::fill(bitmap_.begin(), bitmap_.end(), 0U);
}


auto Bitmap::set(unsigned int const seed_value) -> void
{
  constexpr auto mask_111 = 7U;
  constexpr auto divider = 3U;  // divide by 8
  bitmap_[seed_value >> divider] |= 1U << (seed_value & mask_111);
}


auto Bitmap::collect_set_bits(unsigned int const bound,
                              Span<unsigned int> const destination) const -> unsigned int
{
  assert(destination.size() >= bound);
  assert(bitmap_.size() * 8U >= bound);  // the bitmap has to cover the range asked for
  constexpr auto bits_per_word = 64U;
  constexpr auto bytes_per_word = 8U;

  auto const * const bytes = bitmap_.data();
  auto const whole_words = bound / bits_per_word;
  auto found = 0U;

  for (auto word_number = 0U; word_number < whole_words; ++word_number)
    {
      /* std::memcpy, not a cast of the byte pointer: reading eight unsigned
         chars as one uint64_t through a punned pointer is what the aliasing
         rules forbid, and the buffer carries no alignment guarantee either. It
         compiles to a single load. */
      uint64_t word = 0;
      std::memcpy(&word,
                  std::next(bytes, static_cast<std::ptrdiff_t>(bytes_per_word) * word_number),
                  sizeof word);
      auto const first_bit = bits_per_word * word_number;
      while (word != 0)
        {
          destination[found] = first_bit + countr_zero(word);
          ++found;
          word &= word - 1;  // clear the lowest set bit
        }
    }

  /* the bits past the last whole word: the buffer is a whole number of bytes,
     not of 64-bit words, so the tail cannot be read the same way */
  for (auto bit = bits_per_word * whole_words; bit < bound; ++bit)
    {
      if (is_set(bit))
        {
          destination[found] = bit;
          ++found;
        }
    }

  return found;
}
