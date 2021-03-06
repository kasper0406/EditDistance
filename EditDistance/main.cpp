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
  
  Test::TestSuite::run_tests();
  
  {
    auto parameter_test_fib = [] (double xfactor, double ffactor) {
      const uint64_t max_seq_len = 610;
      const uint64_t fib_len = 15;
      const int trials = 1;
      
      return [xfactor, ffactor, max_seq_len, fib_len, trials] () -> int {
        Benchmark::run_benchmark<Compression::LCSBlowUpAligner<SLP::LZSLPBuilder, DIST::MergingDISTRepository<DIST::PermutationDISTTable, DIST::PermutationLCSMerger>>>(trials, xfactor, ffactor, Benchmark::FastaInput("hg_repetitive.fa", max_seq_len), Benchmark::FastaInput("hg_repetitive.fa", max_seq_len));
        
        Benchmark::run_benchmark<Compression::LCSBlowUpAligner<SLP::LZSLPBuilder, DIST::MergingDISTRepository<DIST::PermutationDISTTable, DIST::PermutationLCSMerger>>>(trials, xfactor, ffactor, Benchmark::FastaInput("hg_nonrepetitive.fa", max_seq_len), Benchmark::FastaInput("hg_nonrepetitive.fa", max_seq_len));
        
        Benchmark::run_benchmark<Compression::LCSBlowUpAligner<SLP::LZSLPBuilder, DIST::MergingDISTRepository<DIST::PermutationDISTTable, DIST::PermutationLCSMerger>>>(trials, xfactor, ffactor, Benchmark::FastaInput("hg_combined.fa", max_seq_len), Benchmark::FastaInput("hg_combined.fa", max_seq_len));
        
        int result = 0;
        
        Compression::LCSBlowUpAligner<SLP::LZSLPBuilder, DIST::MergingDISTRepository<DIST::PermutationDISTTable, DIST::PermutationLCSMerger>> fibAligner(Benchmark::fib_string(fib_len), Benchmark::fib_string(fib_len), xfactor, ffactor);
        result ^= fibAligner.edit_distance();
        
        Compression::LCSBlowUpAligner<SLP::LZSLPBuilder, DIST::MergingDISTRepository<DIST::PermutationDISTTable, DIST::PermutationLCSMerger>> randomAligner(Benchmark::generate_string(max_seq_len), Benchmark::generate_string(max_seq_len), xfactor, ffactor);
        result ^= randomAligner.edit_distance();
        
        return result;
      };
    };
    
    auto parameter_test_random = [] (double xfactor, double ffactor) {
      const uint64_t max_seq_len = 610;
      const uint64_t fib_len = 15;
      const int trials = 1;
      
      return [xfactor, ffactor, max_seq_len, fib_len, trials] () -> int {
        int result = 0;
        
        Compression::LCSBlowUpAligner<SLP::LZSLPBuilder, DIST::MergingDISTRepository<DIST::PermutationDISTTable, DIST::PermutationLCSMerger>> randomAligner(Benchmark::generate_string(max_seq_len), Benchmark::generate_string(max_seq_len), xfactor, ffactor);
        result ^= randomAligner.edit_distance();
        
        return result;
      };
    };
    
    ofstream out_fib("parameter-test-fib.dat");
    Benchmark::test_partition_constants(parameter_test_fib, out_fib);
    out_fib.close();
    
    ofstream out_random("parameter-test-random.dat");
    Benchmark::test_partition_constants(parameter_test_random, out_random);
    out_random.close();
  }
  
  {
    const uint16_t trials = 5;
    const uint64_t max_len = 4000000;
    
    Benchmark::benchmark_compression(trials, Benchmark::FibonacciInput(30));
    Benchmark::benchmark_compression(trials, Benchmark::UniformRandomInput(max_len));
    
    Benchmark::benchmark_compression(trials, Benchmark::FastaInput("hg_repetitive.fa", max_len));
    Benchmark::benchmark_compression(trials, Benchmark::FastaInput("hg_nonrepetitive.fa", max_len));
    Benchmark::benchmark_compression(trials, Benchmark::FastaInput("hg_combined.fa", max_len));
    
    Benchmark::benchmark_min_multiply(trials);
    Benchmark::benchmark_max_multiply(trials);
    Benchmark::benchmark_slow_max_multiply(trials);
    
    const uint64_t max_seq_len = 13000;
    const uint64_t fib_len = 30;
    
    // const double xfactor = 0.7;
    // const double ffactor = 2.5;
    
    const double xfactor_fib = 2.0;
    const double ffactor_fib = 2.0;
    
    const double xfactor_random = 0.4;
    const double ffactor_random = 2.0;
    
    cout << "HG repetitive" << endl;
    Benchmark::run_benchmark<Simple::EditDistance>(trials, xfactor_random, ffactor_random, Benchmark::FastaInput("hg_repetitive.fa", max_seq_len), Benchmark::FastaInput("hg_repetitive.fa", max_seq_len));
    Benchmark::run_benchmark<Compression::LCSBlowUpAligner<SLP::LZSLPBuilder, DIST::MergingDISTRepository<DIST::PermutationDISTTable, DIST::PermutationLCSMerger>>>
    (trials, xfactor_random, ffactor_random, Benchmark::FastaInput("hg_repetitive.fa", max_seq_len), Benchmark::FastaInput("hg_repetitive.fa", max_seq_len));
    
    cout << "HG nonrepetitive" << endl;
    Benchmark::run_benchmark<Simple::EditDistance>(trials, xfactor_random, ffactor_random, Benchmark::FastaInput("hg_nonrepetitive.fa", max_seq_len), Benchmark::FastaInput("hg_nonrepetitive.fa", max_seq_len));
    Benchmark::run_benchmark<Compression::LCSBlowUpAligner<SLP::LZSLPBuilder, DIST::MergingDISTRepository<DIST::PermutationDISTTable, DIST::PermutationLCSMerger>>>
    (trials, xfactor_random, ffactor_random, Benchmark::FastaInput("hg_nonrepetitive.fa", max_seq_len), Benchmark::FastaInput("hg_nonrepetitive.fa", max_seq_len));
    
    cout << "HG combined" << endl;
    Benchmark::run_benchmark<Simple::EditDistance>(trials, xfactor_random, ffactor_random, Benchmark::FastaInput("hg_combined.fa", max_seq_len), Benchmark::FastaInput("hg_combined.fa", max_seq_len));
    Benchmark::run_benchmark<Compression::LCSBlowUpAligner<SLP::LZSLPBuilder, DIST::MergingDISTRepository<DIST::PermutationDISTTable, DIST::PermutationLCSMerger>>>
    (trials, xfactor_random, ffactor_random, Benchmark::FastaInput("hg_combined.fa", max_seq_len), Benchmark::FastaInput("hg_combined.fa", max_seq_len));

    cout << "Fibonacci" << endl;
    Benchmark::run_benchmark<Compression::LCSBlowUpAligner<SLP::LZSLPBuilder, DIST::MergingDISTRepository<DIST::PermutationDISTTable, DIST::PermutationLCSMerger>>>
    (trials, xfactor_fib, ffactor_fib, Benchmark::FibonacciInput(fib_len), Benchmark::FibonacciInput(fib_len));
    Benchmark::run_benchmark<Simple::EditDistance>(trials, xfactor_fib, ffactor_fib, Benchmark::FibonacciInput(fib_len - 6), Benchmark::FibonacciInput(fib_len - 6));

    cout << "Random strings" << endl;
    Benchmark::run_benchmark<Simple::EditDistance>(trials, xfactor_random, ffactor_random, Benchmark::UniformRandomInput(max_seq_len), Benchmark::UniformRandomInput(max_seq_len));
    Benchmark::run_benchmark<Compression::LCSBlowUpAligner<SLP::LZSLPBuilder, DIST::MergingDISTRepository<DIST::PermutationDISTTable, DIST::PermutationLCSMerger>>>
      (trials, xfactor_random, ffactor_random, Benchmark::UniformRandomInput(max_seq_len), Benchmark::UniformRandomInput(max_seq_len));
  }

#ifdef USE_COUNTING_ALLOCATOR
  cout << "Max allocated: " << CountingAllocator::max_allocated << endl;
  cout << "Currently allocated: " << CountingAllocator::currently_allocated << endl;
#endif
  
  return 0;
}
