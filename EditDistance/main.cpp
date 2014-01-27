#include <iostream>
#include <vector>
#include <memory>

#include "utils/allocator.hpp"
#ifdef USE_COUNTING_ALLOCATOR
  template <class T> using StlAllocator = GenericStlAllocator<T, CountingAllocator>;
#else
  template <class T> using StlAllocator = std::allocator<T>;
#endif

#include "simple/simple.hpp"
#include "benchmark/benchmarker.hpp"

using namespace std;

int main(int argc, char* argv[])
{
  {
    uint16_t trials = 5;
    Benchmark::run_benchmark<Simple>(trials);
  }

#ifdef USE_COUNTING_ALLOCATOR
  cout << "Max allocated: " << CountingAllocator::max_allocated << endl;
  cout << "Currently allocated: " << CountingAllocator::currently_allocated << endl;
#endif
  
  return 0;
}
