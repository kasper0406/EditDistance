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
  
  /*
  {
    const string A = "abbababbbabab";
    const string B = "babaaababaab";
    const int64_t x = 1;
    
    unique_ptr<SLP::SLP> slpA = SLP::SimpleCompressionSLPBuilder::build(A);
    unique_ptr<SLP::SLP> slpB = SLP::SimpleCompressionSLPBuilder::build(B);
    auto partitionA = SLP::Partitioner::partition(*slpA, x);
    auto partitionB = SLP::Partitioner::partition(*slpB, x);
    
    DIST::LCSDISTRepository<DIST::SimpleLCSDISTTable> distRepo(*slpA, *slpB);
    distRepo.build(get<1>(partitionA), get<1>(partitionB));
    
    for (int64_t i = 0; i < get<1>(partitionA).size(); ++i) {
      for (int64_t j = 0; j < get<1>(partitionB).size(); ++j) {
        cout << distRepo(get<1>(partitionA)[i], get<1>(partitionB)[j]).matrix() << endl;
        
        if (get<1>(partitionA)[i]->associatedString.length() == 1 &&
            get<1>(partitionB)[j]->associatedString.length() == 1) {
          cout << "Both length 1: " << (get<1>(partitionA)[i]->associatedString == get<1>(partitionB)[j]->associatedString) << endl;
        }
      }
    }
  }
   */
  
  {
    const int64_t x = 20;
    const uint16_t trials = 1;
    Benchmark::run_benchmark<Simple::EditDistance>(trials, x);
    Benchmark::run_benchmark<Compression::EditDistanceAligner<SLP::SimpleCompressionSLPBuilder, DIST::SimpleDISTRepository<DIST::EditDistanceDISTTable>>>(trials, x);
    Benchmark::run_benchmark<Compression::EditDistanceAligner<SLP::SimpleCompressionSLPBuilder, DIST::MergingDISTRepository<DIST::EditDistanceDISTTable, DIST::SimpleMerger>>>(trials, x);
    Benchmark::run_benchmark<Compression::LCSBlowUpAligner<SLP::SimpleCompressionSLPBuilder, DIST::LCSDISTRepository<DIST::SimpleLCSDISTTable>>>(trials, x);
  }

#ifdef USE_COUNTING_ALLOCATOR
  cout << "Max allocated: " << CountingAllocator::max_allocated << endl;
  cout << "Currently allocated: " << CountingAllocator::currently_allocated << endl;
#endif
  
  return 0;
}
