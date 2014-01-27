#pragma once

#include <cstdint>
#include <cassert>

class CountingAllocator {
public:
  static int64_t allocated;
  static int64_t currently_allocated;
  static int64_t max_allocated;
  
  static int64_t local_peek;
  
  template<typename T> static T* allocate(size_t n) {
    allocated += sizeof(T) * n;
    currently_allocated += sizeof(T) * n;
    if (max_allocated < currently_allocated) {
      max_allocated = currently_allocated;
    }
    if (local_peek < currently_allocated) {
      local_peek = currently_allocated;
    }
    
    return new T[n];
  }
  
  template<typename T> static void deallocate(T* ptr, size_t n) {
    currently_allocated -= sizeof(T) * n;
    
    delete[] ptr;
  }
  
  static void reset_local_peek() {
    local_peek = 0;
  }
};

int64_t CountingAllocator::allocated = 0;
int64_t CountingAllocator::currently_allocated = 0;
int64_t CountingAllocator::max_allocated = 0;
int64_t CountingAllocator::local_peek = 0;

template<typename T, class Alloc>
class GenericStlAllocator {
public:
  typedef T        value_type;
  typedef T*       pointer;
  typedef const T* const_pointer;
  typedef T&       reference;
  typedef const T& const_reference;
  typedef std::size_t    size_type;
  typedef std::ptrdiff_t difference_type;
  
  pointer allocate(size_t n) {
    return Alloc::template allocate<T>(n);
  }
  
  void deallocate(pointer ptr, size_t n) {
    Alloc::template deallocate<T>(ptr, n);
  }
};