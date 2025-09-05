// Copyright (c) 2011 Google, Inc.
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
//
// CityHash, by Geoff Pike and Jyrki Alakuijala
//
// This file provides CityHash64() and related functions.
//
// It's probably possible to create even faster hash functions by
// writing a program that systematically explores some of the space of
// possible hash functions, by using SIMD instructions, or by
// compromising on hash quality.

#include "config.h"
#include <city.h>
#include "utils/os_byteswap.hpp"

#include <algorithm>  // std::swap
#include <cstdint>  // int32_t, uint8_t, uint32_t, uint64_t
#include <cstdio>  // std::size_t
#include <cstring>  // std::memcpy, std::memset
#include <utility> // std::pair, std::make_pair


// anonymous namespace: limit visibility and usage to this translation unit
namespace {

  // Some primes between 2^63 and 2^64 for various uses.
  constexpr uint64_t k0 = 0xc3a5c85c97cb3127ULL;  // uint128 seeds?
  constexpr uint64_t k1 = 0xb492b66fbe98f273ULL;  // uint128 seeds?
  constexpr uint64_t k2 = 0x9ae16a3b2f90404fULL;  // uint64 initialization value?


  inline auto Uint128Low64(const uint128 & a_pair) -> uint64_t {
    return a_pair.first;
  }


  inline auto Uint128High64(const uint128 & a_pair) -> uint64_t {
    return a_pair.second;
  }


  auto unaligned_load64(const char * src) -> uint64_t {
    uint64_t result = 0;
    std::memcpy(&result, src, sizeof(result));
    return result;
  }


  auto unaligned_load32(const char * src) -> uint32_t {
    uint32_t result = 0;
    std::memcpy(&result, src, sizeof(result));
    return result;
  }


  // C++23 refactoring:
  //   #include <bit>
  //   constexpr auto uint32_in_expected_order(uint32_t x) noexcept -> uint32_t {
  //     if constexpr (std::endian::native == std::endian::big) {
  //         return std::byteswap(x); }
  //     else {
  //         return x; }
  // }
  constexpr auto uint32_in_expected_order(uint32_t src) noexcept -> uint32_t {
#ifdef WORDS_BIGENDIAN
    return bswap_32(src);
#else
    return src;
#endif
  }


  constexpr auto uint64_in_expected_order(uint64_t src) noexcept -> uint64_t {
#ifdef WORDS_BIGENDIAN
    return bswap_64(src);
#else
    return src;
#endif
  }


  auto Fetch64(const char * src) -> uint64_t {
    return uint64_in_expected_order(unaligned_load64(src));
  }


  auto Fetch32(const char * src) -> uint32_t {
    return uint32_in_expected_order(unaligned_load32(src));
  }


  // Bitwise right rotate.  Normally this will compile to a single
  // instruction, especially if the shift is a manifest constant.
  // C++20 refactoring: std::rotr()
  auto Rotate(uint64_t val, int shift) -> uint64_t {
    static constexpr auto sixtyfour = 64U;
    // Avoid shifting by 64: doing so yields an undefined result.
    return shift == 0 ? val : ((val >> shift) | (val << (sixtyfour - shift)));
  }


  auto ShiftMix(uint64_t val) -> uint64_t {
    static constexpr auto fourtyseven = 47U;
    return val ^ (val >> fourtyseven);
  }


  // Hash 128 input bits down to 64 bits of output.
  // This is intended to be a reasonably good hash function.
  inline auto Hash128to64(const uint128 & a_pair) -> uint64_t {
    // Murmur-inspired hashing.
    static constexpr auto divider = 47U;
    static constexpr uint64_t kMul = 0x9ddfea08eb382d69ULL;
    uint64_t lower_part = (Uint128Low64(a_pair) ^ Uint128High64(a_pair)) * kMul;
    lower_part ^= (lower_part >> divider);
    uint64_t higher_part = (Uint128High64(a_pair) ^ lower_part) * kMul;
    higher_part ^= (higher_part >> divider);
    higher_part *= kMul;
    return higher_part;
  }


  auto HashLen16(uint64_t u, uint64_t v) -> uint64_t {
    return Hash128to64(uint128(u, v));
  }


  auto HashLen16(uint64_t u, uint64_t v, uint64_t mul) -> uint64_t {
    // Murmur-inspired hashing.
    uint64_t a = (u ^ v) * mul;
    a ^= (a >> 47U);
    uint64_t b = (v ^ a) * mul;
    b ^= (b >> 47U);
    b *= mul;
    return b;
  }


  auto HashLen0to16(const char * seq, std::size_t len) -> uint64_t {
    if (len >= 8) {
      const uint64_t mul = k2 + (len * 2);
      const uint64_t a = Fetch64(seq) + k2;
      const uint64_t b = Fetch64(seq + len - 8);
      const uint64_t c = (Rotate(b, 37) * mul) + a;
      const uint64_t d = (Rotate(a, 25) + b) * mul;
      return HashLen16(c, d, mul);
    }
    if (len >= 4) {
      const uint64_t mul = k2 + (len * 2);
      const uint64_t a = Fetch32(seq);
      return HashLen16(len + (a << 3U), Fetch32(seq + len - 4), mul);
    }
    if (len > 0) {
      uint8_t const a = seq[0];
      uint8_t const b = seq[len >> 1U];
      uint8_t const c = seq[len - 1];
      const uint32_t y = static_cast<uint32_t>(a) + (static_cast<uint32_t>(b) << 8U);
      const uint32_t z = len + (static_cast<uint32_t>(c) << 2U);
      return ShiftMix((y * k2) ^ (z * k0)) * k2;
    }
    return k2;  // initialization value for empty sequences
  }


  // This probably works well for 16-byte strings as well, but it may be overkill
  // in that case.
  auto HashLen17to32(const char * seq, std::size_t len) -> uint64_t {
    const uint64_t mul = k2 + (len * 2);
    const uint64_t a = Fetch64(seq) * k1;
    const uint64_t b = Fetch64(seq + 8);
    const uint64_t c = Fetch64(seq + len - 8) * mul;
    const uint64_t d = Fetch64(seq + len - 16) * k2;
    return HashLen16(Rotate(a + b, 43) + Rotate(c, 30) + d,
                     a + Rotate(b + k2, 18) + c, mul);
  }


  // Return a 16-byte hash for 48 bytes.  Quick and dirty.
  // Callers do best to use "random-looking" values for a and b.
  auto WeakHashLen32WithSeeds(uint64_t w, uint64_t x, uint64_t y, uint64_t z,
                              uint64_t a, uint64_t b)
    -> std::pair<uint64_t, uint64_t> {
    a += w;
    b = Rotate(b + a + z, 21);
    const uint64_t c = a;
    a += x;
    a += y;
    b += Rotate(a, 44);
    return std::make_pair(a + z, b + c);
  }


  // Return a 16-byte hash for s[0] ... s[31], a, and b.  Quick and dirty.
  auto WeakHashLen32WithSeeds(const char * seq, uint64_t a, uint64_t b)
    -> std::pair<uint64_t, uint64_t> {
    return WeakHashLen32WithSeeds(Fetch64(seq),
                                  Fetch64(seq + 8),
                                  Fetch64(seq + 16),
                                  Fetch64(seq + 24),
                                  a,
                                  b);
  }


  // Return an 8-byte hash for 33 to 64 bytes.
  auto HashLen33to64(const char * seq, std::size_t len) -> uint64_t {
    const uint64_t mul = k2 + (len * 2);
    uint64_t a = Fetch64(seq) * k2;
    uint64_t b = Fetch64(seq + 8);
    const uint64_t c = Fetch64(seq + len - 24);
    const uint64_t d = Fetch64(seq + len - 32);
    const uint64_t e = Fetch64(seq + 16) * k2;
    const uint64_t f = Fetch64(seq + 24) * 9;
    const uint64_t g = Fetch64(seq + len - 8);
    const uint64_t h = Fetch64(seq + len - 16) * mul;
    const uint64_t u = Rotate(a + g, 43) + ((Rotate(b, 30) + c) * 9);
    const uint64_t v = ((a + g) ^ d) + f + 1;
    const uint64_t w = bswap_64((u + v) * mul) + h;
    const uint64_t x = Rotate(e + f, 42) + c;
    const uint64_t y = (bswap_64((v + w) * mul) + g) * mul;
    const uint64_t z = e + f + c;
    a = bswap_64(((x + z) * mul) + y) + b;
    b = ShiftMix(((z + a) * mul) + d + h) * mul;
    return b + x;
  }


  // A subroutine for CityHash128().  Returns a decent 128-bit hash for strings
  // of any length representable in signed long.  Based on City and Murmur.
  auto CityMurmur(const char * seq, std::size_t len, uint128 seed) -> uint128 {
    uint64_t a = Uint128Low64(seed);
    uint64_t b = Uint128High64(seed);
    uint64_t c = 0;
    uint64_t d = 0;
    signed long l = len - 16;
    if (l <= 0) {  // len <= 16
      a = ShiftMix(a * k1) * k1;
      c = (b * k1) + HashLen0to16(seq, len);
      d = ShiftMix(a + (len >= 8 ? Fetch64(seq) : c));
    } else {  // len > 16
      c = HashLen16(Fetch64(seq + len - 8) + k1, a);
      d = HashLen16(b + len, c + Fetch64(seq + len - 16));
      a += d;
      do {
        a ^= ShiftMix(Fetch64(seq) * k1) * k1;
        a *= k1;
        b ^= a;
        c ^= ShiftMix(Fetch64(seq + 8) * k1) * k1;
        c *= k1;
        d ^= c;
        seq += 16;
        l -= 16;
      } while (l > 0);
    }
    a = HashLen16(a, c);
    b = HashLen16(d, b);
    return uint128(a ^ b, HashLen16(b, a));
  }


  constexpr auto likely(bool condition) noexcept -> bool {
#if HAVE_BUILTIN_EXPECT
    static constexpr auto is_expected_to_be_true = 1L;
    // !!condition converts to a strict boolean value
    return __builtin_expect(!!condition, is_expected_to_be_true);
#else
    return condition;
#endif
  }


  auto CityHash128WithSeed(const char * seq, std::size_t len, uint128 seed) -> uint128 {
    if (len < 128) {
      return CityMurmur(seq, len, seed);
    }

    // We expect len >= 128 to be the common case.  Keep 56 bytes of state:
    // v, w, x, y, and z.
    std::pair<uint64_t, uint64_t> v;
    std::pair<uint64_t, uint64_t> w;
    uint64_t x = Uint128Low64(seed);
    uint64_t y = Uint128High64(seed);
    uint64_t z = len * k1;
    v.first = (Rotate(y ^ k1, 49) * k1) + Fetch64(seq);
    v.second = (Rotate(v.first, 42) * k1) + Fetch64(seq + 8);
    w.first = (Rotate(y + z, 35) * k1) + x;
    w.second = Rotate(x + Fetch64(seq + 88), 53) * k1;

    // This is the same inner loop as CityHash64(), manually unrolled.
    do {
      x = Rotate(x + y + v.first + Fetch64(seq + 8), 37) * k1;
      y = Rotate(y + v.second + Fetch64(seq + 48), 42) * k1;
      x ^= w.second;
      y += v.first + Fetch64(seq + 40);
      z = Rotate(z + w.first, 33) * k1;
      v = WeakHashLen32WithSeeds(seq, v.second * k1, x + w.first);
      w = WeakHashLen32WithSeeds(seq + 32, z + w.second, y + Fetch64(seq + 16));
      std::swap(z, x);
      seq += 64;
      x = Rotate(x + y + v.first + Fetch64(seq + 8), 37) * k1;
      y = Rotate(y + v.second + Fetch64(seq + 48), 42) * k1;
      x ^= w.second;
      y += v.first + Fetch64(seq + 40);
      z = Rotate(z + w.first, 33) * k1;
      v = WeakHashLen32WithSeeds(seq, v.second * k1, x + w.first);
      w = WeakHashLen32WithSeeds(seq + 32, z + w.second, y + Fetch64(seq + 16));
      std::swap(z, x);
      seq += 64;
      len -= 128;
    } while (likely(len >= 128));  // hot path
    x += Rotate(v.first + z, 49) * k0;
    y = (y * k0) + Rotate(w.second, 37);
    z = (z * k0) + Rotate(w.first, 27);
    w.first *= 9;
    v.first *= k0;
    // If 0 < len < 128, hash up to 4 chunks of 32 bytes each from the end of s.
    for (std::size_t tail_done = 0; tail_done < len; ) {
      tail_done += 32;
      y = (Rotate(x + y, 42) * k0) + v.second;
      w.first += Fetch64(seq + len - tail_done + 16);
      x = (x * k0) + w.first;
      z += w.second + Fetch64(seq + len - tail_done);
      w.second += v.first;
      v = WeakHashLen32WithSeeds(seq + len - tail_done, v.first + z, v.second);
      v.first *= k0;
    }
    // At this point our 56 bytes of state should contain more than
    // enough information for a strong 128-bit hash.  We use two
    // different 56-byte-to-8-byte hashes to get a 16-byte final result.
    x = HashLen16(x, v.first);
    y = HashLen16(y + z, w.first);
    return uint128(HashLen16(x + v.second, w.second) + y,
                   HashLen16(x + w.second, y + v.second));
  }

}  // end of anonymous namespace


auto CityHash64(const char * seq, std::size_t len) -> uint64_t {
  if (len <= 16) {
    return HashLen0to16(seq, len);
  }
  if (len <= 32) {
    return HashLen17to32(seq, len);
  }
  if (len <= 64) {
    return HashLen33to64(seq, len);
  }

  // For strings over 64 bytes we hash the end first, and then as we
  // loop we keep 56 bytes of state: v, w, x, y, and z.
  uint64_t x = Fetch64(seq + len - 40);
  uint64_t y = Fetch64(seq + len - 16) + Fetch64(seq + len - 56);
  uint64_t z = HashLen16(Fetch64(seq + len - 48) + len, Fetch64(seq + len - 24));
  auto v = WeakHashLen32WithSeeds(seq + len - 64, len, z);
  auto w = WeakHashLen32WithSeeds(seq + len - 32, y + k1, x);
  x = (x * k1) + Fetch64(seq);

  // Decrease len to the nearest multiple of 64, and operate on 64-byte chunks.
  len = (len - 1) & ~static_cast<std::size_t>(63);
  do {
    x = Rotate(x + y + v.first + Fetch64(seq + 8), 37) * k1;
    y = Rotate(y + v.second + Fetch64(seq + 48), 42) * k1;
    x ^= w.second;
    y += v.first + Fetch64(seq + 40);
    z = Rotate(z + w.first, 33) * k1;
    v = WeakHashLen32WithSeeds(seq, v.second * k1, x + w.first);
    w = WeakHashLen32WithSeeds(seq + 32, z + w.second, y + Fetch64(seq + 16));
    std::swap(z, x);
    seq += 64;
    len -= 64;
  } while (len != 0);
  return HashLen16(HashLen16(v.first, w.first) + (ShiftMix(y) * k1) + z,
                   HashLen16(v.second, w.second) + x);
}


auto CityHash128(const char * seq, std::size_t len) -> uint128 {
  return len >= 16 ?
      CityHash128WithSeed(seq + 16, len - 16,
                          uint128(Fetch64(seq), Fetch64(seq + 8) + k0)) :
      CityHash128WithSeed(seq, len, uint128(k0, k1));
}
