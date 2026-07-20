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

#pragma once

#include "os/system.hpp"  // xmalloc, xfree
#include <cstddef>  // std::size_t


/* A minimal standard-library allocator that obtains memory through xmalloc and
   releases it through xfree. xmalloc reports an out-of-memory condition through
   fatal() (a clean message, then exit) and never returns null, so an exhausted
   allocation ends the program exactly as the raw db buffers always did, instead
   of throwing std::bad_alloc -- which, under -fno-exceptions, would
   std::terminate. Stateless: any two instances compare equal, so a container
   may move or swap its storage freely. C++11 std::allocator_traits synthesises
   the rest of the allocator interface (rebind, construct, max_size, ...). */
template <class Element>
struct FatalAllocator
{
  using value_type = Element;

  FatalAllocator() = default;

  // enables a container to rebind the allocator for its internal element type
  template <class Other>
  FatalAllocator(FatalAllocator<Other> const & /*other*/) noexcept {}  // NOLINT: implicit by design

  auto allocate(std::size_t const count) -> Element *
  {
    // a container never requests more than max_size() elements, so the product
    // cannot overflow std::size_t; xmalloc fatals rather than returning null
    return static_cast<Element *>(xmalloc(count * sizeof(Element)));
  }

  auto deallocate(Element * const pointer, std::size_t const /*count*/) noexcept -> void
  {
    xfree(static_cast<void *>(pointer));
  }
};


template <class Left, class Right>
auto operator==(FatalAllocator<Left> const & /*lhs*/, FatalAllocator<Right> const & /*rhs*/) noexcept -> bool
{
  return true;
}

template <class Left, class Right>
auto operator!=(FatalAllocator<Left> const & /*lhs*/, FatalAllocator<Right> const & /*rhs*/) noexcept -> bool
{
  return false;
}
