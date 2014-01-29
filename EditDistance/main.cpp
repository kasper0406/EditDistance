#include <iostream>
#include <vector>
#include <memory>
#include <string>
#include <map>
#include <array>

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

using namespace std;

/**
 * NOTICE: Strings are NOT generated uniformly at random!!!
 */
string generate_string(uint64_t length, vector<char> alphabet = { 'a', 'c', 'g', 't' })
{
  stringstream ss;
  
  for (uint64_t i = 0; i < length; ++i)
    ss << alphabet[rand() % alphabet.size()]; // Not fair randomness, but its good enough for this testing!
  
  return ss.str();
}

void test_partition_generation() {
  uint64_t x = 5;
  
  for (uint64_t length = 1; length < 1000; length = max((uint64_t)((double)length * 1.2), length + 1)) {
    string input = generate_string(length);
    
    for (uint64_t x = 1; x < input.size(); x = max((uint64_t)((double)x * 1.2), x + 1)) {
      SLP::SLP slp = SLP::SimpleCompressionSLPBuilder::build(input);
      // SLP::SLP slp = SLP::SimpleSLPBuilder::build(input);
      
      /*
      if (x == 1) {
        cout << SLP::SLPUnfoldedPrinter::toDot(slp) << endl;
        cout << SLP::SLPFoldedPrinter::toDot(slp) << endl;
      }
       */
      
      auto partition = SLP::BlockConstructor::buildBlocks(slp, x);
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
    uint16_t trials = 5;
    Benchmark::run_benchmark<Simple>(trials);
  }
  
  {
    test_partition_generation();
  }

#ifdef USE_COUNTING_ALLOCATOR
  cout << "Max allocated: " << CountingAllocator::max_allocated << endl;
  cout << "Currently allocated: " << CountingAllocator::currently_allocated << endl;
#endif
  
  return 0;
}
