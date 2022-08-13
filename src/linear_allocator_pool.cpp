
/*
Adapted by Rob Egan from https://github.com/mtrebi/memory-allocators

MIT License

Copyright (c) 2016 Mariano Trebino

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include "linear_allocator_pool.hpp"

#include <stdlib.h> /* malloc, free */
#include <cassert>

#include "upcxx_utils/log.hpp"

using upcxx_utils::get_size_str;

LinearAllocatorPool::LinearAllocatorPool(const std::size_t block_size, const std::size_t total_size)
    : m_block_size(block_size)
    , m_totalSize{0}
    , m_used{0}
    , m_peak{0}
    , m_num_allocations{0} {
  Reserve(total_size);
  assert(m_totalSize >= total_size);
}

LinearAllocatorPool::~LinearAllocatorPool() { FreeBlocks(); }

const std::size_t LinearAllocatorPool::CalculatePadding(const std::size_t baseAddress, const std::size_t alignment) {
  const std::size_t multiplier = (baseAddress / alignment) + 1;
  const std::size_t alignedAddress = multiplier * alignment;
  const std::size_t padding = alignedAddress - baseAddress;
  return padding;
}

void *LinearAllocatorPool::Allocate(const std::size_t size, const std::size_t alignment) {
  std::size_t padding = 0;
  std::size_t paddedAddress = 0;

  if (size > m_block_size) {
    // Allocate its own entire block at the front and return that
    LOG("Allocated more than blocksize=", m_block_size, " size=", size, "\n");
    auto m_start_ptr = malloc(size);
    if (m_start_ptr == nullptr) throw std::bad_alloc();
    m_blocks.push_front({m_start_ptr, size});
    m_used += size;
    m_peak = std::max(m_peak, m_used);
    m_num_allocations++;
    return m_start_ptr;
  }

  if (m_blocks.empty() || m_blocks.back().second + size > m_block_size) {
    // if the previous one is unused, move this filled block in the back up to the front
    if (m_blocks.size() > 1 && m_blocks[m_blocks.size() - 2].second == 0) {
      m_blocks.push_front(m_blocks.back());
      m_blocks.pop_back();
    } else {
      AddBlock();
    }
    assert(!m_blocks.empty());
    assert(m_blocks.back().second == 0);
  }

  auto m_start_ptr = m_blocks.back().first;
  auto m_offset = m_blocks.back().second;

  const std::size_t currentAddress = (std::size_t)m_start_ptr + m_offset;

  if (alignment != 0 && m_offset % alignment != 0) {
    // Alignment is required. Find the next aligned memory address and update offset
    padding = CalculatePadding(currentAddress, alignment);
    if (m_offset + padding + size > m_block_size) {
      AddBlock();
      m_start_ptr = m_blocks.back().first;
      m_offset = m_blocks.back().second;
      assert(m_offset == 0);
      padding = 0;
    }
  }

  m_offset += padding;
  const std::size_t nextAddress = currentAddress + padding;
  m_offset += size;

  m_used += size + padding;
  m_peak = std::max(m_peak, m_used);
  m_num_allocations++;

  m_blocks.back().second = m_offset;

  return (void *)nextAddress;
}

void LinearAllocatorPool::Free(void *ptr) { assert(false && "Use Reset() method, not Free"); }

void LinearAllocatorPool::Init() {
  FreeBlocks();
  auto sz = m_blocks.size() * m_block_size;
  while (sz < m_totalSize) {
    AddBlock();
    sz += m_block_size;
  }
}

void LinearAllocatorPool::AddBlock() {
  auto m_start_ptr = calloc(1,m_block_size); // allocate and touch
  if (m_start_ptr == nullptr) throw std::bad_alloc();
  m_blocks.push_back({m_start_ptr, 0});
  LOG("Added new ", get_size_str(m_block_size), " block. size=", m_blocks.size(), " num_allocations=", m_num_allocations,
      " m_used=", get_size_str(m_used), " m_peak=", get_size_str(m_peak), " allocator=", (size_t)this, "\n");
}

void LinearAllocatorPool::Reset() {
  FreeBlocks();
  m_used = 0;
  m_peak = 0;
}

void LinearAllocatorPool::FreeBlocks() {
  bool freed = false;
  for (auto &block : m_blocks) {
    freed = true;
    free(block.first);
    block.first = nullptr;
  }
  if (freed) {
    LOG("Freed blocks size=", m_blocks.size(), " num_allocations=", m_num_allocations, " m_used=", get_size_str(m_used),
        " m_peak=", get_size_str(m_peak), " allocator=", (size_t)this, "\n");
    Blocks().swap(m_blocks);
    assert(m_blocks.empty());
  }
}

void LinearAllocatorPool::Reserve(const std::size_t size) {
  if (size > m_totalSize) {
    const std::size_t n_blocks = (size + m_block_size - 1) / m_block_size;
    m_totalSize = n_blocks * m_block_size;
    Init();
  }
}