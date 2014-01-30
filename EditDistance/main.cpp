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

void test_partition_generation() {
  for (uint64_t length = 1; length < 1000; length = max((uint64_t)((double)length * 1.2), length + 1)) {
    string input = Benchmark::generate_string(length);
    
    for (uint64_t x = 1; x < input.size(); x = max((uint64_t)((double)x * 1.2), x + 1)) {
      SLP::SLP slp = SLP::SimpleCompressionSLPBuilder::build(input);
      // SLP::SLP slp = SLP::SimpleSLPBuilder::build(input);
      
      /*
      if (x == 1) {
        cout << SLP::SLPUnfoldedPrinter::toDot(slp) << endl;
        cout << SLP::SLPFoldedPrinter::toDot(slp) << endl;
      }
       */
      
      auto partition = get<0>(SLP::Partitioner::partition(slp, x));
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
  /*
  {
    Simple simple("abaab", "aaa");
    cout << simple.edit_distance() << endl;
    
    uint16_t trials = 5;
    Benchmark::run_benchmark<Simple>(trials);
  }
   */
  
  {
    test_partition_generation();
    
    const string A = Benchmark::generate_string(1000);
    const string B = Benchmark::generate_string(1000);
    
    /*
    const string A = "abbabaabababaaabbbbbbaaababbbbababbbbbabbbabbbb";
    const string B = "ababaaaabbbbababaaababbbbbababbbaaaba";
     */
    
    {
      Simple simple(A, B);
      cout << "Should be: " << simple.edit_distance() << endl;
    }
    
    const uint64_t x = 15;
    
    SLP::SLP slpA = SLP::SimpleCompressionSLPBuilder::build(A);
    vector<SLP::Production*> partitionA, blocksA;
    tie(partitionA, blocksA) = SLP::Partitioner::partition(slpA, x);
    
    SLP::SLP slpB = SLP::SimpleCompressionSLPBuilder::build(B);
    vector<SLP::Production*> partitionB, blocksB;
    tie(partitionB, blocksB) = SLP::Partitioner::partition(slpB, x);
    
    DIST::SimpleDISTRepository DISTRepo(blocksA, blocksB);
    
    // Base case for column
    vector<uint64_t> prev_column(A.size() + 1, 0), cur_column(A.size() + 1, 0);
    for (uint64_t Apos = 0; Apos <= A.size(); ++Apos) {
      prev_column[Apos] = Apos;
    }
    
    int64_t Bpos = 0;
    for (auto b : partitionB) {
      const int64_t b_len = b->associatedString.size();
      
      vector<uint64_t> block_row(b_len + 1, 0);
      for (int64_t i = 0; i <= b_len; i++)
        block_row[i] = Bpos + i;
      
      int64_t Apos = 0;
      for (auto a : partitionA) {
        const int64_t a_len = a->associatedString.size();
        assert(block_row[0] == prev_column[Apos]);
        
        // Build inputs to block
        vector<uint64_t> I(a_len + b_len + 1, 0);
        for (int64_t i = 0; i <= a_len; ++i) {
          I[i] = prev_column[Apos + a_len - i];
          // cout << "I[" << i << "] = " << I[i] << endl;
        }
        for (int64_t i = a_len + 1; i < I.size(); ++i)
          I[i] = block_row[i - a_len];
        
        // Computes outputs from block
        auto O = DISTRepo.apply(a, b, I);
        /*
        cout << "I: " << "\t";
        for (int64_t i = 0; i < I.size(); ++i)
          cout << I[i] << "\t";
        cout << endl;
        
        cout << "O: " << "\t";
        for (int64_t i = 0; i < O.size(); ++i)
          cout << O[i] << "\t";
        cout << endl;
         */
        
        assert(O.size() == I.size());
        // Write output
        for (int64_t i = 0; i <= b_len; ++i)
          block_row[i] = O[i];
        for (int64_t i = b_len; i < O.size(); ++i)
          cur_column[Apos + O.size() - i - 1] = O[i];
        
        Apos += a_len;
      }
      
      swap(prev_column, cur_column);
      
      Bpos += b_len;
      if (Bpos == B.size())
        cout << "Found dist: " << block_row[block_row.size() - 1] << endl;
    }
  }

#ifdef USE_COUNTING_ALLOCATOR
  cout << "Max allocated: " << CountingAllocator::max_allocated << endl;
  cout << "Currently allocated: " << CountingAllocator::currently_allocated << endl;
#endif
  
  return 0;
}
