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

#include <algorithm>  // std::swap
#include <cstdint>  // int32_t, uint8_t, uint32_t, uint64_t
#include <cstdio>  // std::size_t
#include <cstring>  // std::memcpy, std::memset
#include <utility> // std::pair, std::make_pair


static auto UNALIGNED_LOAD64(const char *p) -> uint64 {
  uint64 result = 0;
  std::memcpy(&result, p, sizeof(result));
  return result;
}

static auto UNALIGNED_LOAD32(const char *p) -> uint32 {
  uint32 result = 0;
  std::memcpy(&result, p, sizeof(result));
  return result;
}

#ifdef _MSC_VER

#include <stdlib.h>
#define bswap_32(x) _byteswap_ulong(x)
#define bswap_64(x) _byteswap_uint64(x)

#elif defined(__APPLE__)

// Mac OS X / Darwin features
#include <libkern/OSByteOrder.h>
#define bswap_32(x) OSSwapInt32(x)
#define bswap_64(x) OSSwapInt64(x)

#elif defined(__FreeBSD__)

#include <sys/endian.h>
#define bswap_32(x) bswap32(x)
#define bswap_64(x) bswap64(x)

#elif defined(__NetBSD__)

#include <sys/types.h>
#include <machine/bswap.h>
#if defined(__BSWAP_RENAME) && !defined(__bswap_32)
#define bswap_32(x) bswap32(x)
#define bswap_64(x) bswap64(x)
#endif

#else

#include <byteswap.h>

#endif

#ifdef WORDS_BIGENDIAN
#define uint32_in_expected_order(x) (bswap_32(x))
#define uint64_in_expected_order(x) (bswap_64(x))
#else
#define uint32_in_expected_order(x) (x)
#define uint64_in_expected_order(x) (x)
#endif

#if !defined(LIKELY)
#if HAVE_BUILTIN_EXPECT
#define LIKELY(x) (__builtin_expect(!!(x), 1))
#else
#define LIKELY(x) (x)
#endif
#endif

static auto Fetch64(const char *p) -> uint64 {
  return uint64_in_expected_order(UNALIGNED_LOAD64(p));
}

static auto Fetch32(const char *p) -> uint32 {
  return uint32_in_expected_order(UNALIGNED_LOAD32(p));
}

// Some primes between 2^63 and 2^64 for various uses.
constexpr uint64 k0 = 0xc3a5c85c97cb3127ULL;
constexpr uint64 k1 = 0xb492b66fbe98f273ULL;
constexpr uint64 k2 = 0x9ae16a3b2f90404fULL;

// Magic numbers for 32-bit hashing.  Copied from Murmur3.
constexpr uint32_t c1 = 0xcc9e2d51;
constexpr uint32_t c2 = 0x1b873593;

// A 32-bit to 32-bit integer hash copied from Murmur3.
static auto fmix(uint32 hash_value) -> uint32
{
  hash_value ^= hash_value >> 16U;
  hash_value *= 0x85ebca6b;
  hash_value ^= hash_value >> 13U;
  hash_value *= 0xc2b2ae35;
  hash_value ^= hash_value >> 16U;
  return hash_value;
}

static auto Rotate32(uint32 val, int shift) -> uint32 {
  // Avoid shifting by 32: doing so yields an undefined result.
  return shift == 0 ? val : ((val >> shift) | (val << (32 - shift)));
}

#undef PERMUTE3
#define PERMUTE3(a, b, c) do { std::swap(a, b); std::swap(a, c); } while (0)

static auto Mur(uint32 a_value, uint32 hash_value) -> uint32 {
  // Helper from Murmur3 for combining two 32-bit values.
  a_value *= c1;
  a_value = Rotate32(a_value, 17);
  a_value *= c2;
  hash_value ^= a_value;
  hash_value = Rotate32(hash_value, 19);
  return (hash_value * 5) + 0xe6546b64;
}

static auto Hash32Len13to24(const char * seq, std::size_t len) -> uint32 {
  const uint32 a = Fetch32(seq - 4 + (len >> 1U));
  const uint32 b = Fetch32(seq + 4);
  const uint32 c = Fetch32(seq + len - 8);
  const uint32 d = Fetch32(seq + (len >> 1U));
  const uint32 e = Fetch32(seq);
  const uint32 f = Fetch32(seq + len - 4);
  const uint32 h = len;

  return fmix(Mur(f, Mur(e, Mur(d, Mur(c, Mur(b, Mur(a, h)))))));
}

static auto Hash32Len0to4(const char * seq, std::size_t len) -> uint32 {
  uint32 b = 0;
  uint32 c = 9;
  for (int i = 0; i < len; i++) {
    const signed char v = seq[i];
    b = b * c1 + v;
    c ^= b;
  }
  return fmix(Mur(b, Mur(len, c)));
}

static auto Hash32Len5to12(const char * seq, std::size_t len) -> uint32 {
  uint32 a = len;
  uint32 b = len * 5;
  uint32 c = 9;
  uint32 const d = b;
  a += Fetch32(seq);
  b += Fetch32(seq + len - 4);
  c += Fetch32(seq + ((len >> 1U) & 4U));
  return fmix(Mur(c, Mur(b, Mur(a, d))));
}

auto CityHash32(const char * seq, std::size_t len) -> uint32 {
  if (len <= 24) {
    return len <= 12 ?
        (len <= 4 ? Hash32Len0to4(seq, len) : Hash32Len5to12(seq, len)) :
        Hash32Len13to24(seq, len);
  }

  // len > 24
  uint32 h = len;
  uint32 g = c1 * len;
  uint32 f = g;
  const uint32 a0 = Rotate32(Fetch32(seq + len - 4) * c1, 17) * c2;
  const uint32 a1 = Rotate32(Fetch32(seq + len - 8) * c1, 17) * c2;
  const uint32 a2 = Rotate32(Fetch32(seq + len - 16) * c1, 17) * c2;
  const uint32 a3 = Rotate32(Fetch32(seq + len - 12) * c1, 17) * c2;
  const uint32 a4 = Rotate32(Fetch32(seq + len - 20) * c1, 17) * c2;
  h ^= a0;
  h = Rotate32(h, 19);
  h = h * 5 + 0xe6546b64;
  h ^= a2;
  h = Rotate32(h, 19);
  h = h * 5 + 0xe6546b64;
  g ^= a1;
  g = Rotate32(g, 19);
  g = g * 5 + 0xe6546b64;
  g ^= a3;
  g = Rotate32(g, 19);
  g = g * 5 + 0xe6546b64;
  f += a4;
  f = Rotate32(f, 19);
  f = f * 5 + 0xe6546b64;
  std::size_t iters = (len - 1) / 20;
  do {
    const uint32 a0 = Rotate32(Fetch32(seq) * c1, 17) * c2;
    const uint32 a1 = Fetch32(seq + 4);
    const uint32 a2 = Rotate32(Fetch32(seq + 8) * c1, 17) * c2;
    const uint32 a3 = Rotate32(Fetch32(seq + 12) * c1, 17) * c2;
    const uint32 a4 = Fetch32(seq + 16);
    h ^= a0;
    h = Rotate32(h, 18);
    h = h * 5 + 0xe6546b64;
    f += a1;
    f = Rotate32(f, 19);
    f = f * c1;
    g += a2;
    g = Rotate32(g, 18);
    g = g * 5 + 0xe6546b64;
    h ^= a3 + a1;
    h = Rotate32(h, 19);
    h = h * 5 + 0xe6546b64;
    g ^= a4;
    g = bswap_32(g) * 5;
    h += a4 * 5;
    h = bswap_32(h);
    f += a0;
    PERMUTE3(f, h, g);
    seq += 20;
  } while (--iters != 0);
  g = Rotate32(g, 11) * c1;
  g = Rotate32(g, 17) * c1;
  f = Rotate32(f, 11) * c1;
  f = Rotate32(f, 17) * c1;
  h = Rotate32(h + g, 19);
  h = h * 5 + 0xe6546b64;
  h = Rotate32(h, 17) * c1;
  h = Rotate32(h + f, 19);
  h = h * 5 + 0xe6546b64;
  h = Rotate32(h, 17) * c1;
  return h;
}

// Bitwise right rotate.  Normally this will compile to a single
// instruction, especially if the shift is a manifest constant.
static auto Rotate(uint64 val, int shift) -> uint64 {
  // Avoid shifting by 64: doing so yields an undefined result.
  return shift == 0 ? val : ((val >> shift) | (val << (64U - shift)));
}

static auto ShiftMix(uint64 val) -> uint64 {
  return val ^ (val >> 47U);
}

static auto HashLen16(uint64 u, uint64 v) -> uint64 {
  return Hash128to64(uint128(u, v));
}

static auto HashLen16(uint64 u, uint64 v, uint64 mul) -> uint64 {
  // Murmur-inspired hashing.
  uint64 a = (u ^ v) * mul;
  a ^= (a >> 47U);
  uint64 b = (v ^ a) * mul;
  b ^= (b >> 47U);
  b *= mul;
  return b;
}

static auto HashLen0to16(const char * seq, std::size_t len) -> uint64 {
  if (len >= 8) {
    const uint64 mul = k2 + (len * 2);
    const uint64 a = Fetch64(seq) + k2;
    const uint64 b = Fetch64(seq + len - 8);
    const uint64 c = (Rotate(b, 37) * mul) + a;
    const uint64 d = (Rotate(a, 25) + b) * mul;
    return HashLen16(c, d, mul);
  }
  if (len >= 4) {
    const uint64 mul = k2 + (len * 2);
    const uint64 a = Fetch32(seq);
    return HashLen16(len + (a << 3U), Fetch32(seq + len - 4), mul);
  }
  if (len > 0) {
    uint8_t const a = seq[0];
    uint8_t const b = seq[len >> 1U];
    uint8_t const c = seq[len - 1];
    const uint32 y = static_cast<uint32>(a) + (static_cast<uint32>(b) << 8U);
    const uint32 z = len + (static_cast<uint32>(c) << 2U);
    return ShiftMix((y * k2) ^ (z * k0)) * k2;
  }
  return k2;
}

// This probably works well for 16-byte strings as well, but it may be overkill
// in that case.
static auto HashLen17to32(const char * seq, std::size_t len) -> uint64 {
  const uint64 mul = k2 + (len * 2);
  const uint64 a = Fetch64(seq) * k1;
  const uint64 b = Fetch64(seq + 8);
  const uint64 c = Fetch64(seq + len - 8) * mul;
  const uint64 d = Fetch64(seq + len - 16) * k2;
  return HashLen16(Rotate(a + b, 43) + Rotate(c, 30) + d,
                   a + Rotate(b + k2, 18) + c, mul);
}

// Return a 16-byte hash for 48 bytes.  Quick and dirty.
// Callers do best to use "random-looking" values for a and b.
static auto WeakHashLen32WithSeeds(uint64 w, uint64 x, uint64 y, uint64 z,
                                   uint64 a, uint64 b)
  -> std::pair<uint64, uint64> {
  a += w;
  b = Rotate(b + a + z, 21);
  const uint64 c = a;
  a += x;
  a += y;
  b += Rotate(a, 44);
  return std::make_pair(a + z, b + c);
}

// Return a 16-byte hash for s[0] ... s[31], a, and b.  Quick and dirty.
static auto WeakHashLen32WithSeeds(const char * seq, uint64 a, uint64 b)
  -> std::pair<uint64, uint64> {
  return WeakHashLen32WithSeeds(Fetch64(seq),
                                Fetch64(seq + 8),
                                Fetch64(seq + 16),
                                Fetch64(seq + 24),
                                a,
                                b);
}

// Return an 8-byte hash for 33 to 64 bytes.
static auto HashLen33to64(const char * seq, std::size_t len) -> uint64 {
  const uint64 mul = k2 + (len * 2);
  uint64 a = Fetch64(seq) * k2;
  uint64 b = Fetch64(seq + 8);
  const uint64 c = Fetch64(seq + len - 24);
  const uint64 d = Fetch64(seq + len - 32);
  const uint64 e = Fetch64(seq + 16) * k2;
  const uint64 f = Fetch64(seq + 24) * 9;
  const uint64 g = Fetch64(seq + len - 8);
  const uint64 h = Fetch64(seq + len - 16) * mul;
  const uint64 u = Rotate(a + g, 43) + ((Rotate(b, 30) + c) * 9);
  const uint64 v = ((a + g) ^ d) + f + 1;
  const uint64 w = bswap_64((u + v) * mul) + h;
  const uint64 x = Rotate(e + f, 42) + c;
  const uint64 y = (bswap_64((v + w) * mul) + g) * mul;
  const uint64 z = e + f + c;
  a = bswap_64(((x + z) * mul) + y) + b;
  b = ShiftMix(((z + a) * mul) + d + h) * mul;
  return b + x;
}

auto CityHash64(const char * seq, std::size_t len) -> uint64 {
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
  uint64 x = Fetch64(seq + len - 40);
  uint64 y = Fetch64(seq + len - 16) + Fetch64(seq + len - 56);
  uint64 z = HashLen16(Fetch64(seq + len - 48) + len, Fetch64(seq + len - 24));
  auto v = WeakHashLen32WithSeeds(seq + len - 64, len, z);
  auto w = WeakHashLen32WithSeeds(seq + len - 32, y + k1, x);
  x = x * k1 + Fetch64(seq);

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

auto CityHash64WithSeed(const char * seq, std::size_t len, uint64 seed) -> uint64 {
  return CityHash64WithSeeds(seq, len, k2, seed);
}

auto CityHash64WithSeeds(const char * seq, std::size_t len,
                         uint64 seed0, uint64 seed1) -> uint64 {
  return HashLen16(CityHash64(seq, len) - seed0, seed1);
}

// A subroutine for CityHash128().  Returns a decent 128-bit hash for strings
// of any length representable in signed long.  Based on City and Murmur.
static auto CityMurmur(const char * seq, std::size_t len, uint128 seed) -> uint128 {
  uint64 a = Uint128Low64(seed);
  uint64 b = Uint128High64(seed);
  uint64 c = 0;
  uint64 d = 0;
  signed long l = len - 16;
  if (l <= 0) {  // len <= 16
    a = ShiftMix(a * k1) * k1;
    c = b * k1 + HashLen0to16(seq, len);
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

auto CityHash128WithSeed(const char * seq, std::size_t len, uint128 seed) -> uint128 {
  if (len < 128) {
    return CityMurmur(seq, len, seed);
  }

  // We expect len >= 128 to be the common case.  Keep 56 bytes of state:
  // v, w, x, y, and z.
  std::pair<uint64, uint64> v;
  std::pair<uint64, uint64> w;
  uint64 x = Uint128Low64(seed);
  uint64 y = Uint128High64(seed);
  uint64 z = len * k1;
  v.first = Rotate(y ^ k1, 49) * k1 + Fetch64(seq);
  v.second = Rotate(v.first, 42) * k1 + Fetch64(seq + 8);
  w.first = Rotate(y + z, 35) * k1 + x;
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
  } while (LIKELY(len >= 128));
  x += Rotate(v.first + z, 49) * k0;
  y = y * k0 + Rotate(w.second, 37);
  z = z * k0 + Rotate(w.first, 27);
  w.first *= 9;
  v.first *= k0;
  // If 0 < len < 128, hash up to 4 chunks of 32 bytes each from the end of s.
  for (std::size_t tail_done = 0; tail_done < len; ) {
    tail_done += 32;
    y = Rotate(x + y, 42) * k0 + v.second;
    w.first += Fetch64(seq + len - tail_done + 16);
    x = x * k0 + w.first;
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

auto CityHash128(const char * seq, std::size_t len) -> uint128 {
  return len >= 16 ?
      CityHash128WithSeed(seq + 16, len - 16,
                          uint128(Fetch64(seq), Fetch64(seq + 8) + k0)) :
      CityHash128WithSeed(seq, len, uint128(k0, k1));
}

#ifdef __SSE4_2__
#include <citycrc.h>
#include <nmmintrin.h>

// Requires len >= 240.
static auto CityHashCrc256Long(const char * s, std::size_t len,
                               uint32 seed, uint64 *result) -> void {
  uint64 a = Fetch64(s + 56) + k0;
  uint64 b = Fetch64(s + 96) + k0;
  uint64 c = result[0] = HashLen16(b, len);
  uint64 d = result[1] = Fetch64(s + 120) * k0 + len;
  uint64 e = Fetch64(s + 184) + seed;
  uint64 f = 0;
  uint64 g = 0;
  uint64 h = c + d;
  uint64 x = seed;
  uint64 y = 0;
  uint64 z = 0;

  // 240 bytes of input per iter.
  std::size_t iters = len / 240;
  len -= iters * 240;
  do {
#undef CHUNK
#define CHUNK(r)                                \
    PERMUTE3(x, z, y);                          \
    b += Fetch64(s);                            \
    c += Fetch64(s + 8);                        \
    d += Fetch64(s + 16);                       \
    e += Fetch64(s + 24);                       \
    f += Fetch64(s + 32);                       \
    a += b;                                     \
    h += f;                                     \
    b += c;                                     \
    f += d;                                     \
    g += e;                                     \
    e += z;                                     \
    g += x;                                     \
    z = _mm_crc32_u64(z, b + g);                \
    y = _mm_crc32_u64(y, e + h);                \
    x = _mm_crc32_u64(x, f + a);                \
    e = Rotate(e, r);                           \
    c += e;                                     \
    s += 40

    CHUNK(0); PERMUTE3(a, h, c);
    CHUNK(33); PERMUTE3(a, h, f);
    CHUNK(0); PERMUTE3(b, h, f);
    CHUNK(42); PERMUTE3(b, h, d);
    CHUNK(0); PERMUTE3(b, h, e);
    CHUNK(33); PERMUTE3(a, h, e);
  } while (--iters > 0);

  while (len >= 40) {
    CHUNK(29);
    e ^= Rotate(a, 20);
    h += Rotate(b, 30);
    g ^= Rotate(c, 40);
    f += Rotate(d, 34);
    PERMUTE3(c, h, g);
    len -= 40;
  }
  if (len > 0) {
    s = s + len - 40;
    CHUNK(33);
    e ^= Rotate(a, 43);
    h += Rotate(b, 42);
    g ^= Rotate(c, 41);
    f += Rotate(d, 40);
  }
  result[0] ^= h;
  result[1] ^= g;
  g += h;
  a = HashLen16(a, g + z);
  x += y << 32;
  b += x;
  c = HashLen16(c, z) + h;
  d = HashLen16(d, e + result[0]);
  g += e;
  h += HashLen16(x, f);
  e = HashLen16(a, d) + g;
  z = HashLen16(b, c) + a;
  y = HashLen16(g, h) + c;
  result[0] = e + z + y + x;
  a = ShiftMix((a + y) * k0) * k0 + b;
  result[1] += a + result[0];
  a = ShiftMix(a * k0) * k0 + c;
  result[2] = a + result[1];
  a = ShiftMix((a + e) * k0) * k0;
  result[3] = a + result[2];
}

// Requires len < 240.
static auto CityHashCrc256Short(const char * s, std::size_t len, uint64 *result) -> void {
  char buf[240];
  std::memcpy(buf, s, len);
  std::memset(buf + len, 0, 240 - len);
  CityHashCrc256Long(buf, 240, ~static_cast<uint32>(len), result);
}

auto CityHashCrc256(const char * s, std::size_t len, uint64 *result) -> void {
  if (LIKELY(len >= 240)) {
    CityHashCrc256Long(s, len, 0, result);
  } else {
    CityHashCrc256Short(s, len, result);
  }
}

auto CityHashCrc128WithSeed(const char * s, std::size_t len, uint128 seed) -> uint128 {
  if (len <= 900) {
    return CityHash128WithSeed(s, len, seed);
  } else {
    uint64 result[4];
    CityHashCrc256(s, len, result);
    uint64 u = Uint128High64(seed) + result[0];
    uint64 v = Uint128Low64(seed) + result[1];
    return uint128(HashLen16(u, v + result[2]),
                   HashLen16(Rotate(v, 32), u * k0 + result[3]));
  }
}

auto CityHashCrc128(const char * s, std::size_t len) -> uint128 {
  if (len <= 900) {
    return CityHash128(s, len);
  } else {
    uint64 result[4];
    CityHashCrc256(s, len, result);
    return uint128(result[2], result[3]);
  }
}

#endif
