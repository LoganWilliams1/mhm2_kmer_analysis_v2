#pragma once

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

#include <cstddef>  // size_t
#include <deque>

class LinearAllocatorPool {
  /*
  Builds a list of fixed size blocks of memory
  When one block fills, a new one is created
  if a larger allocation is requested than the block_size, a new immediately filled block is made for that allocation
  choose the block size wisely as memory can be wasted at the end of the block if the next request does not fit perfectly
  Memory is freed only when this LinearAllocator is deconstructed
  */
 protected:
  std::size_t m_block_size;
  std::size_t m_totalSize;
  std::size_t m_used;
  std::size_t m_peak;
  std::size_t m_num_allocations;
  using Block = std::pair<void *, std::size_t>;
  using Blocks = std::deque<Block>;
  Blocks m_blocks;

 public:
  LinearAllocatorPool(const std::size_t block_size, const std::size_t total_size = 0);

  ~LinearAllocatorPool();

  static const std::size_t CalculatePadding(const std::size_t baseAddress, const std::size_t alignment);

  void *Allocate(const std::size_t size, const std::size_t alignment = 0);

  void Free(void *ptr);

  void Init();

  void AddBlock();

  void Reset();

  void Reserve(const std::size_t size);

 private:
  LinearAllocatorPool(LinearAllocatorPool &linearAllocator) = delete;
  void FreeBlocks();
};
