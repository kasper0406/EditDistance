#include <iostream>
#include <vector>
#include <memory>
#include <string>
#include <map>
#include <array>
#include <memory>

#include "utils/allocator.hpp"
#ifdef USE_COUNTING_ALLOCATOR
  template <class T> using StlAllocator = GenericStlAllocator<T, CountingAllocator>;
#else
  template <class T> using StlAllocator = std::allocator<T>;
#endif

#include "test/tests.hpp"

#include "simple/simple.hpp"
#include "benchmark/benchmarker.hpp"

#include "compression/SLP.hpp"
#include "compression/DIST.hpp"
#include "compression/slpalign.hpp"

using namespace std;
using namespace Compression;

int main(int argc, char* argv[])
{
  Test::TestSuite::run_tests();
  
  {
    const uint64_t x = 20;
    const uint16_t trials = 1;
    Benchmark::run_benchmark<Simple>(trials, x);
    Benchmark::run_benchmark<SLPAlign<SLP::SimpleCompressionSLPBuilder, DIST::SimpleDISTRepository>>(trials, x);
  }

#ifdef USE_COUNTING_ALLOCATOR
  cout << "Max allocated: " << CountingAllocator::max_allocated << endl;
  cout << "Currently allocated: " << CountingAllocator::currently_allocated << endl;
#endif
  
  return 0;
}
