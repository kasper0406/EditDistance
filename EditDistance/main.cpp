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
  
  // cout << DIST::PermutationDISTTable::BaseCase(true)->unfoldH() << endl;
  
  return 0;
  
  {
    const int64_t x = 20;
    const uint16_t trials = 1;
    
    Benchmark::run_benchmark<Compression::LCSBlowUpAligner<SLP::SimpleCompressionSLPBuilder, DIST::MergingDISTRepository<DIST::SimpleLCSDISTTable, DIST::SimpleLCSDISTMerger>>>(trials, x);
    
    /*
    Benchmark::run_benchmark<Simple::EditDistance>(trials, x);
    Benchmark::run_benchmark<Compression::EditDistanceAligner<SLP::SimpleCompressionSLPBuilder, DIST::SimpleDISTRepository<DIST::EditDistanceDISTTable>>>(trials, x);
    Benchmark::run_benchmark<Compression::EditDistanceAligner<SLP::SimpleCompressionSLPBuilder, DIST::MergingDISTRepository<DIST::EditDistanceDISTTable, DIST::SimpleMerger>>>(trials, x);
    Benchmark::run_benchmark<Compression::LCSBlowUpAligner<SLP::SimpleCompressionSLPBuilder, DIST::MergingDISTRepository<DIST::SimpleLCSDISTTable, DIST::SimpleLCSDISTMerger>>>(trials, x);
    Benchmark::run_benchmark<Compression::LCSBlowUpAligner<SLP::SimpleCompressionSLPBuilder, DIST::LCSDISTRepository<DIST::SimpleLCSDISTTable>>>(trials, x);
     */
  }

#ifdef USE_COUNTING_ALLOCATOR
  cout << "Max allocated: " << CountingAllocator::max_allocated << endl;
  cout << "Currently allocated: " << CountingAllocator::currently_allocated << endl;
#endif
  
  return 0;
}
