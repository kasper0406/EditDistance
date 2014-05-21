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
    auto parameter_test = [] (double xfactor, double ffactor) {
      const uint64_t max_seq_len = 610;
      const uint64_t fib_len = 15;
      const int trials = 1;
      
      return [xfactor, ffactor, max_seq_len, fib_len, trials] () -> int {
        /*
        Benchmark::run_benchmark<Compression::LCSBlowUpAligner<SLP::LZSLPBuilder, DIST::MergingDISTRepository<DIST::PermutationDISTTable, DIST::PermutationLCSMerger>>>(trials, xfactor, ffactor, Benchmark::FastaInput("hg_repetitive.fa", max_seq_len), Benchmark::FastaInput("hg_repetitive.fa", max_seq_len));
        
        Benchmark::run_benchmark<Compression::LCSBlowUpAligner<SLP::LZSLPBuilder, DIST::MergingDISTRepository<DIST::PermutationDISTTable, DIST::PermutationLCSMerger>>>(trials, xfactor, ffactor, Benchmark::FastaInput("hg_nonrepetitive.fa", max_seq_len), Benchmark::FastaInput("hg_nonrepetitive.fa", max_seq_len));
        
        Benchmark::run_benchmark<Compression::LCSBlowUpAligner<SLP::LZSLPBuilder, DIST::MergingDISTRepository<DIST::PermutationDISTTable, DIST::PermutationLCSMerger>>>(trials, xfactor, ffactor, Benchmark::FastaInput("hg_combined.fa", max_seq_len), Benchmark::FastaInput("hg_combined.fa", max_seq_len));
         */
        
        int result = 0;
        
        Compression::LCSBlowUpAligner<SLP::LZSLPBuilder, DIST::MergingDISTRepository<DIST::PermutationDISTTable, DIST::PermutationLCSMerger>> fibAligner(Benchmark::fib_string(fib_len), Benchmark::fib_string(fib_len), xfactor, ffactor);
        result ^= fibAligner.edit_distance();
        
        Compression::LCSBlowUpAligner<SLP::LZSLPBuilder, DIST::MergingDISTRepository<DIST::PermutationDISTTable, DIST::PermutationLCSMerger>> randomAligner(Benchmark::generate_string(max_seq_len), Benchmark::generate_string(max_seq_len), xfactor, ffactor);
        result ^= randomAligner.edit_distance();
        
        return result;
      };
    };
    
    ofstream out("parameter-test.dat");
    Benchmark::test_partition_constants(parameter_test, out);
    out.close();
  }
  
  return 0;
  
  {
    const double xfactor = 2;
    const uint16_t trials = 10;
    const uint64_t max_len = 4000000;
    
    /*
    Benchmark::benchmark_compression(trials, Benchmark::FibonacciInput(30));
    Benchmark::benchmark_compression(trials, Benchmark::UniformRandomInput(max_len));
     */
    
    /*
    Benchmark::benchmark_compression(trials, Benchmark::FastaInput("hg_repetitive.fa", max_len));
    Benchmark::benchmark_compression(trials, Benchmark::FastaInput("hg_nonrepetitive.fa", max_len));
    Benchmark::benchmark_compression(trials, Benchmark::FastaInput("hg_combined.fa", max_len));
     */
    
    /*
    Benchmark::benchmark_compression(trials, Benchmark::FastaInput("genome1.fa", max_len));
    Benchmark::benchmark_compression(trials, Benchmark::FastaInput("genome2.fa", max_len));
    Benchmark::benchmark_compression(trials, Benchmark::FibonacciInput(30));
    Benchmark::benchmark_compression(trials, Benchmark::UniformRandomInput(max_len));
     */
    
    /*
    Benchmark::benchmark_min_multiply(trials);
    Benchmark::benchmark_max_multiply(trials);
    Benchmark::benchmark_slow_max_multiply(trials);
     */
    
    Benchmark::benchmark_max_multiply(trials);
    
    /*
    const uint64_t max_seq_len = 13000;
    const uint64_t fib_len = 30;
    
    cout << "HG repetitive" << endl;
    Benchmark::run_benchmark<Simple::EditDistance>(trials, xfactor, Benchmark::FastaInput("hg_repetitive.fa", max_seq_len), Benchmark::FastaInput("hg_repetitive.fa", max_seq_len));
    Benchmark::run_benchmark<Compression::LCSBlowUpAligner<SLP::LZSLPBuilder, DIST::MergingDISTRepository<DIST::PermutationDISTTable, DIST::PermutationLCSMerger>>>
    (trials, xfactor, Benchmark::FastaInput("hg_repetitive.fa", max_seq_len), Benchmark::FastaInput("hg_repetitive.fa", max_seq_len));
    
    cout << "HG nonrepetitive" << endl;
    Benchmark::run_benchmark<Simple::EditDistance>(trials, xfactor, Benchmark::FastaInput("hg_nonrepetitive.fa", max_seq_len), Benchmark::FastaInput("hg_nonrepetitive.fa", max_seq_len));
    Benchmark::run_benchmark<Compression::LCSBlowUpAligner<SLP::LZSLPBuilder, DIST::MergingDISTRepository<DIST::PermutationDISTTable, DIST::PermutationLCSMerger>>>
    (trials, xfactor, Benchmark::FastaInput("hg_nonrepetitive.fa", max_seq_len), Benchmark::FastaInput("hg_nonrepetitive.fa", max_seq_len));
    
    cout << "HG combined" << endl;
    Benchmark::run_benchmark<Simple::EditDistance>(trials, xfactor, Benchmark::FastaInput("hg_combined.fa", max_seq_len), Benchmark::FastaInput("hg_combined.fa", max_seq_len));
    Benchmark::run_benchmark<Compression::LCSBlowUpAligner<SLP::LZSLPBuilder, DIST::MergingDISTRepository<DIST::PermutationDISTTable, DIST::PermutationLCSMerger>>>
    (trials, xfactor, Benchmark::FastaInput("hg_combined.fa", max_seq_len), Benchmark::FastaInput("hg_combined.fa", max_seq_len));
     */
    
    /*
    cout << "Fibonacci" << endl;
    Benchmark::run_benchmark<Simple::EditDistance>(trials, xfactor, Benchmark::FibonacciInput(fib_len - 6), Benchmark::FibonacciInput(fib_len - 6));
    Benchmark::run_benchmark<Compression::LCSBlowUpAligner<SLP::LZSLPBuilder, DIST::MergingDISTRepository<DIST::PermutationDISTTable, DIST::PermutationLCSMerger>>>
      (trials, xfactor, Benchmark::FibonacciInput(fib_len), Benchmark::FibonacciInput(fib_len));

    cout << "Random strings" << endl;
    Benchmark::run_benchmark<Simple::EditDistance>(trials, xfactor, Benchmark::UniformRandomInput(max_seq_len), Benchmark::UniformRandomInput(max_seq_len));
    Benchmark::run_benchmark<Compression::LCSBlowUpAligner<SLP::LZSLPBuilder, DIST::MergingDISTRepository<DIST::PermutationDISTTable, DIST::PermutationLCSMerger>>>
      (trials, xfactor, Benchmark::UniformRandomInput(max_seq_len), Benchmark::UniformRandomInput(max_seq_len));
     */
    
     
     
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
