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

#include "core/minheap.hpp"
#include <algorithm>  // std::sort
#include <cstddef>  // std::size_t


/* implement a priority queue with a min heap binary array structure */
/* elements with the lowest count should be at the top (root) */

// refactoring note: std::priority_queue is a poor fit here. It offers no
// bounded (top-k) insertion, so keeping only the best 'capacity' elements
// would cost an extra sift per candidate, and it hides its container, so the
// sort-then-drain-best-first access below (sort() + pop_last()) is not
// expressible. A small hand-written heap keeps both cheap.

/*
  To keep track of the n best potential target sequences, we store
  them in a min heap. The root element corresponds to the least good
  target, while the best elements are found at the leaf nodes. This
  makes it simple to decide whether a new target should be included or
  not, because it just needs to be compared to the root note.  The
  list will be fully sorted before use when we want to find the best
  element and then the second best and so on.
*/

namespace {

auto elem_less(elem_t const & lhs, elem_t const & rhs) -> bool
{
  /* first: lower count, larger length, lower seqno */
  if (lhs.count != rhs.count)
    {
      return lhs.count < rhs.count;
    }
  if (lhs.length != rhs.length)
    {
      return lhs.length > rhs.length;
    }
  return lhs.seqno > rhs.seqno;
}

}  // namespace


Minheap::Minheap(int capacity)
  : capacity_(static_cast<std::size_t>(capacity))
{
  array_.reserve(capacity_);
}


auto Minheap::is_empty() const -> bool
{
  return array_.empty();
}


auto Minheap::clear() -> void
{
  array_.clear();
}


auto Minheap::replace_root(elem_t tmp) -> void
{
  /* remove the element at the root, then swap children up
     to the root and insert tmp at suitable place */

  /* start with root */
  std::size_t parent = 0;
  std::size_t nth_child = (2 * parent) + 1;
  auto const count = array_.size();

  /* while at least one child */
  while (nth_child < count)
    {
      /* if two children: swap with the one with smallest value */
      if ((nth_child + 1 < count) and
          elem_less(array_[nth_child + 1], array_[nth_child]))
        {
          ++nth_child;
        }

      /* swap parent and child if child has lower value */
      if (elem_less(array_[nth_child], tmp))
        {
          array_[parent] = array_[nth_child];
        }
      else
        {
          break;
        }

      /* step down */
      parent = nth_child;
      nth_child = (2 * parent) + 1;
    }

  array_[parent] = tmp;
}


auto Minheap::add(elem_t const & element) -> void
{
  if (array_.size() < capacity_)
    {
      /* space for another item at end; swap upwards */

      array_.push_back(element);
      auto index = array_.size() - 1;
      while (index > 0)
        {
          auto const pos = (index - 1) / 2;
          if (not elem_less(element, array_[pos]))
            {
              break;
            }
          array_[index] = array_[pos];
          index = pos;
        }
      array_[index] = element;
    }
  else if ((not array_.empty()) and elem_less(array_.front(), element))
    {
      /* replace the root if new element is larger than root */
      replace_root(element);
    }
}


auto Minheap::sort() -> void
{
  std::sort(array_.begin(), array_.end(), elem_less);
}


auto Minheap::pop_last() -> elem_t
{
  /* return top element and restore order */
  if (array_.empty())
    {
      return elem_t{0, 0, 0};
    }
  elem_t const last = array_.back();
  array_.pop_back();
  return last;
}
