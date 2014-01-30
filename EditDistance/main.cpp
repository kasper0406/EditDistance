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

#include "simple/simple.hpp"
#include "benchmark/benchmarker.hpp"

#include "compression/SLP.hpp"
#include "compression/DIST.hpp"
#include "compression/slpalign.hpp"

using namespace std;
using namespace Compression;

void test_partition_generation() {
  for (uint64_t length = 1; length < 1000; length = max((uint64_t)((double)length * 1.2), length + 1)) {
    string input = Benchmark::generate_string(length);
    
    for (uint64_t x = 1; x < input.size(); x = max((uint64_t)((double)x * 1.2), x + 1)) {
      unique_ptr<SLP::SLP> slp = SLP::SimpleCompressionSLPBuilder::build(input);
      // SLP::SLP slp = SLP::SimpleSLPBuilder::build(input);
      
      /*
      if (x == 1) {
        cout << SLP::SLPUnfoldedPrinter::toDot(slp) << endl;
        cout << SLP::SLPFoldedPrinter::toDot(slp) << endl;
      }
       */
      
      auto partition = get<0>(SLP::Partitioner::partition(*slp, x));
      stringstream partitionString;
      for (auto block : partition) {
        // cout << "X" << dynamic_cast<SLP::NonTerminal*>(block)->name() << "\t" << block->associatedString << endl;
        partitionString << block->associatedString;
      }
      
      // cout << length << "\t" << x << "\t" << partitionString.str() << endl;
      if (partitionString.str() != input)
        throw runtime_error("Partition string different from input!");
    }
  }
}

int main(int argc, char* argv[])
{
  {
    const uint64_t x = 20;
    const uint16_t trials = 5;
    Benchmark::run_benchmark<Simple>(trials, x);
    Benchmark::run_benchmark<SLPAlign<SLP::SimpleCompressionSLPBuilder, DIST::SimpleDISTRepository>>(trials, x);
  }
  
  /*
  {
    test_partition_generation();
    
    const string A = Benchmark::generate_string(1000);
    const string B = Benchmark::generate_string(1000);
    
//    const string A = "abbabaabababaaabbbbbbaaababbbbababbbbbabbbabbbb";
//    const string B = "ababaaaabbbbababaaababbbbbababbbaaaba";
    
    {
      const uint64_t x = 5;
      SLPAlign<SLP::SimpleCompressionSLPBuilder, DIST::SimpleDISTRepository> aligner(A, B, x);
      cout << "Found: " << aligner.edit_distance() << endl;
    }
    
    {
      Simple simple(A, B);
      cout << "Should be: " << simple.edit_distance() << endl;
    }
  }
  */

#ifdef USE_COUNTING_ALLOCATOR
  cout << "Max allocated: " << CountingAllocator::max_allocated << endl;
  cout << "Currently allocated: " << CountingAllocator::currently_allocated << endl;
#endif
  
  return 0;
}
