#include <iostream>
#include <vector>
#include <memory>
#include <string>
#include <map>
#include <array>
#include <memory>

#include <sdsl/suffix_trees.hpp>

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
using namespace sdsl;

int main(int argc, char* argv[])
{
#ifndef NDEBUG
  cout << "Running in DEBUG MODE!" << endl << endl;
#else
  cout << "Running in RELEASE MODE!" << endl << endl;
#endif
  
  // TODO: Write test for blow up slp
  Test::TestSuite::run_tests();
  
  {
    const double xfactor = 3;
    const uint16_t trials = 1;
    
    /*
    Benchmark::benchmark_compression(trials, Benchmark::FastaInput("genome1.fa", max_len));
    Benchmark::benchmark_compression(trials, Benchmark::FastaInput("genome2.fa", max_len));
    Benchmark::benchmark_compression(trials, Benchmark::FibonacciInput(30));
    Benchmark::benchmark_compression(trials, Benchmark::UniformRandomInput(max_len));
     */
    
    // Benchmark::benchmark_min_multiply(trials);
    
    // Benchmark::benchmark_max_multiply(trials);
    // Benchmark::benchmark_slow_max_multiply(trials);
    
    const uint64_t max_seq_len = 50000;
    const uint64_t fib_len = 30;
    
    cout << "DNA sequences" << endl;
    Benchmark::run_benchmark<Simple::EditDistance>(trials, xfactor, Benchmark::FastaInput("genome1.fa", max_seq_len), Benchmark::FastaInput("genome2.fa", max_seq_len));
    Benchmark::run_benchmark<Compression::LCSBlowUpAligner<SLP::LZSLPBuilder, DIST::MergingDISTRepository<DIST::PermutationDISTTable, DIST::PermutationLCSMerger>>>
      (trials, xfactor, Benchmark::FastaInput("genome1.fa", max_seq_len), Benchmark::FastaInput("genome2.fa", max_seq_len));
     
    cout << "Fibonacci" << endl;
    Benchmark::run_benchmark<Simple::EditDistance>(trials, xfactor, Benchmark::FibonacciInput(fib_len - 10), Benchmark::FibonacciInput(fib_len - 10));
    Benchmark::run_benchmark<Compression::LCSBlowUpAligner<SLP::LZSLPBuilder, DIST::MergingDISTRepository<DIST::PermutationDISTTable, DIST::PermutationLCSMerger>>>
      (trials, xfactor, Benchmark::FibonacciInput(fib_len), Benchmark::FibonacciInput(fib_len));

    cout << "Random strings" << endl;
    Benchmark::run_benchmark<Simple::EditDistance>(trials, xfactor, Benchmark::UniformRandomInput(max_seq_len), Benchmark::UniformRandomInput(max_seq_len));
    Benchmark::run_benchmark<Compression::LCSBlowUpAligner<SLP::LZSLPBuilder, DIST::MergingDISTRepository<DIST::PermutationDISTTable, DIST::PermutationLCSMerger>>>
      (trials, xfactor, Benchmark::UniformRandomInput(max_seq_len), Benchmark::UniformRandomInput(max_seq_len));
    
    // Benchmark::run_benchmark<Compression::LCSBlowUpAligner<SLP::SimpleCompressionSLPBuilder, DIST::MergingDISTRepository<DIST::SimpleLCSDISTTable, DIST::SimpleLCSDISTMerger>>>(trials, x);
    
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
