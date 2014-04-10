#pragma once

#include <vector>
#include <functional>
#include <memory>
#include <utility>
#include <string>

#include "../simple/simple.hpp"

#include "../compression/SLP.hpp"
#include "../compression/DIST.hpp"
#include "../compression/maxmultiply.hpp"
#include "../compression/slpalign.hpp"

#include "../benchmark/benchmarker.hpp"

#include "../utils/unionfind.hpp"

namespace Test {
  using namespace std;
  using namespace Compression;
  
  class TestSuite {
  public:
    static void run_tests() {
      const int64_t colwidth = 80;
      
      vector<pair<string, function<bool()>>> tests = {
        { "Test LZ factorization", test_lz_factorization },
        
        { "Testing union find", test_union_find },
        { "Testing interval union find", test_interval_union_find },
        
        { "Verify result of sample using LCS blow up SLP and merging DIST.", test_slp_compression_align<LCSBlowUpAligner<SLP::SimpleCompressionSLPBuilder, DIST::MergingDISTRepository<DIST::PermutationDISTTable, DIST::PermutationLCSMerger>>> },
        
        { "Verify result of sample using LCS blow up LZ-SLP and merging DIST.", test_slp_compression_align<LCSBlowUpAligner<SLP::LZSLPBuilder, DIST::MergingDISTRepository<DIST::PermutationDISTTable, DIST::PermutationLCSMerger>>> },
        
        { "Partition generation with trivial SLP.", test_partition_generation<SLP::SimpleSLPBuilder> },
        { "Partition generation with simple compression SLP.", test_partition_generation<SLP::SimpleCompressionSLPBuilder> },
        { "Test simple horizontal DIST merge.", test_horizontal_merge<DIST::SimpleDISTRepository<DIST::EditDistanceDISTTable>, DIST::SimpleMerger> },
        { "Test simple vertical DIST merge.", test_vertical_merge<DIST::SimpleDISTRepository<DIST::EditDistanceDISTTable>, DIST::SimpleMerger> },
        { "Test construction of DIST tables using simple merging.", test_dist_repository<DIST::SimpleDISTRepository<DIST::EditDistanceDISTTable>, DIST::MergingDISTRepository<DIST::EditDistanceDISTTable, DIST::SimpleMerger>, SLP::SimpleCompressionSLPBuilder> },
        { "Verify result of sample using simple compression SLP and simple DIST.", test_slp_compression_align<EditDistanceAligner<SLP::SimpleCompressionSLPBuilder, DIST::SimpleDISTRepository<DIST::EditDistanceDISTTable>>> },
        { "Verify result of sample using simple compression SLP and simple merging DIST.", test_slp_compression_align<EditDistanceAligner<SLP::SimpleCompressionSLPBuilder, DIST::MergingDISTRepository<DIST::EditDistanceDISTTable, DIST::SimpleMerger>>> },
        { "Test blow-up method for computing edit-distance from LCS.", test_blow_up_method },
        { "Test LCS vertical simple DIST merge.", test_vertical_merge<DIST::LCSDISTRepository<DIST::SimpleLCSDISTTable>, DIST::SimpleLCSDISTMerger> },
        { "Test LCS horizontal simple DIST merge.", test_horizontal_merge<DIST::LCSDISTRepository<DIST::SimpleLCSDISTTable>, DIST::SimpleLCSDISTMerger> },
        { "Test construction of DIST tables using simple merging LCS DIST", test_dist_repository<DIST::LCSDISTRepository<DIST::SimpleLCSDISTTable>, DIST::MergingDISTRepository<DIST::SimpleLCSDISTTable, DIST::SimpleLCSDISTMerger>> },
        { "Verify result of sample using LCS blow up SLP and simple DIST.", test_slp_compression_align<LCSBlowUpAligner<SLP::SimpleCompressionSLPBuilder, DIST::LCSDISTRepository<DIST::SimpleLCSDISTTable>>> },        
        { "Verify result of sample using LCS blow up SLP and simple merging DIST.", test_slp_compression_align<LCSBlowUpAligner<SLP::SimpleCompressionSLPBuilder, DIST::MergingDISTRepository<DIST::SimpleLCSDISTTable, DIST::SimpleLCSDISTMerger>>> },
        { "Test merging permutation LCS DIST builder.", test_dist_repository<DIST::LCSDISTRepository<DIST::SimpleLCSDISTTable>, DIST::MergingDISTRepository<DIST::PermutationDISTTable, DIST::PermutationLCSMerger>> },
        
        { "Partition generation with LZ SLP.", test_partition_generation<SLP::LZSLPBuilder> },
        { "Test LZ SLP builder", test_slp_builder<SLP::LZSLPBuilder> }
      };
      
      bool success = true;
      cout << "Running test suite..." << endl;
      for (auto test : tests) {
        cout << left << setw(colwidth) << get<0>(test) << flush;
        if (get<1>(test)()) {
          cout << "PASSED" << endl;
        } else {
          cout << "FAILED!" << endl;
          success = false;
        }
      }
      
      cout << endl;
      if (success)
        cout << "All test-cases passed :)!" << endl;
      else
        cout << "Oh no :(! There a test case is not working." << endl;
      cout << endl << endl;
    }
    
  private:
    template <class SLPBuilder>
    static bool test_slp_builder() {
      for (int64_t length = 10; length < 1000; length = max((int64_t)((double)length * 1.2), length + 1)) {
        string input = Benchmark::generate_string(length);
        auto tree = SLP::LZSLPBuilder::build(input);
        
        assert(tree->derivedLength() == input.length());
        if (SLP::StringDeriver::getDerivedString(tree->root()) != input)
          return false;
      }
      
      return true;
    }
    
    template <class SLPBuilder>
    static bool test_partition_generation() {
      for (int64_t length = 10; length < 1000; length = max((int64_t)((double)length * 1.4), length + 1)) {
        string input = Benchmark::generate_string(length, { 'a', 'b', 'c' });
        
        for (int64_t x = 5; x < input.size(); x = max((int64_t)((double)x * 4), x + 1)) {
          unique_ptr<SLP::SLP> slp = SLPBuilder::build(input);
          
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
          if (partitionString.str() != input) {
            cout << "Was: " << partitionString.str() << endl;
            cout << "Should be: " << input << endl;
            
            return false;
          }
        }
      }
      
      return true;
    }
    
    template <class DISTRepo, class Merger>
    static bool test_horizontal_merge()
    {
       /*
      const string A = "abbaccabcbbcbabcbacbccacbabcbacacccacbbbbaabaac";
      const string B = "aacabbacbbcc";
      const string C = "cabc";
       */
      
      const uint64_t N = 8;
      
      for (int64_t i = 1; i < N; ++i) {
        for (int64_t j = 1; j < N; ++j) {
          for (int64_t k = 1; k < N; ++k) {
            const string A = Benchmark::generate_string(i);
            const string B = Benchmark::generate_string(j);
            const string C = Benchmark::generate_string(k);
            
            unique_ptr<SLP::SLP> slpA = SLP::SimpleCompressionSLPBuilder::build(A);
            unique_ptr<SLP::SLP> slpB = SLP::SimpleCompressionSLPBuilder::build(B);
            auto partitionA = SLP::Partitioner::partition(*slpA, A.size());
            auto partitionB = SLP::Partitioner::partition(*slpB, B.size());
            
            DISTRepo distRepoAB(*slpA, *slpB);
            distRepoAB.build(get<1>(partitionA), get<1>(partitionB));
            
            unique_ptr<SLP::SLP> slpA_ = SLP::SimpleCompressionSLPBuilder::build(A);
            unique_ptr<SLP::SLP> slpC = SLP::SimpleCompressionSLPBuilder::build(C);
            auto partitionA_ = SLP::Partitioner::partition(*slpA_, A.size());
            auto partitionC = SLP::Partitioner::partition(*slpC, C.size());
            
            DISTRepo distRepoAC(*slpA_, *slpC);
            distRepoAC.build(get<1>(partitionA_), get<1>(partitionC));
            
            //  cout << get<2>(distRepoAB(slpA->root(), slpB->root())) << endl;
            //  cout << get<2>(distRepoAC(slpA_->root(), slpC->root())) << endl;
            
            unique_ptr<SLP::SLP> slpA__ = SLP::SimpleCompressionSLPBuilder::build(A);
            unique_ptr<SLP::SLP> slpBC = SLP::SimpleCompressionSLPBuilder::build(B + C);
            auto partitionA__ = SLP::Partitioner::partition(*slpA__, A.size());
            auto partitionBC = SLP::Partitioner::partition(*slpBC, (B + C).size());
            
            DISTRepo distRepoABC(*slpA__, *slpBC);
            distRepoABC.build(get<1>(partitionA__), get<1>(partitionBC));
            
            auto merged = Merger::horizontalMerge(&distRepoAB(slpA->root(), slpB->root()),
                                                  &distRepoAC(slpA_->root(), slpC->root()));
      
            if (distRepoABC(slpA__->root(), slpBC->root()) != *merged)
              return false;
          }
        }
      }
      
      return true;
    }
    
    template <class DISTRepo, class Merger>
    static bool test_vertical_merge()
    {
      /*
      const string A = "aacabbacbbcccabc";
      const string B = "abbaccabcbbcbabcbacb";
      const string C = "ccacbabcbacacccacbbbbaabaac";
       */
      
      const uint64_t N = 8;
      
      for (int64_t i = 1; i < N; ++i) {
        for (int64_t j = 1; j < N; ++j) {
          for (int64_t k = 1; k < N; ++k) {
            const string A = Benchmark::generate_string(i);
            const string B = Benchmark::generate_string(j);
            const string C = Benchmark::generate_string(k);
            
            unique_ptr<SLP::SLP> slpA = SLP::SimpleCompressionSLPBuilder::build(A);
            unique_ptr<SLP::SLP> slpB = SLP::SimpleCompressionSLPBuilder::build(B);
            auto partitionA = SLP::Partitioner::partition(*slpA, A.size());
            auto partitionB = SLP::Partitioner::partition(*slpB, B.size());
            
            DISTRepo distRepoAB(*slpB, *slpA);
            distRepoAB.build(get<1>(partitionB), get<1>(partitionA));
            
            unique_ptr<SLP::SLP> slpA_ = SLP::SimpleCompressionSLPBuilder::build(A);
            unique_ptr<SLP::SLP> slpC = SLP::SimpleCompressionSLPBuilder::build(C);
            auto partitionA_ = SLP::Partitioner::partition(*slpA_, A.size());
            auto partitionC = SLP::Partitioner::partition(*slpC, C.size());
            
            DISTRepo distRepoAC(*slpC, *slpA_);
            distRepoAC.build(get<1>(partitionC), get<1>(partitionA_));
            
            //  cout << get<2>(distRepoAB(slpB->root(), slpA->root())) << endl;
            //  cout << get<2>(distRepoAC(slpC->root(), slpA_->root())) << endl;
            
            unique_ptr<SLP::SLP> slpA__ = SLP::SimpleCompressionSLPBuilder::build(A);
            unique_ptr<SLP::SLP> slpBC = SLP::SimpleCompressionSLPBuilder::build(B + C);
            auto partitionA__ = SLP::Partitioner::partition(*slpA__, A.size());
            auto partitionBC = SLP::Partitioner::partition(*slpBC, (B + C).size());
            
            DISTRepo distRepoABC(*slpBC, *slpA__);
            distRepoABC.build(get<1>(partitionBC), get<1>(partitionA__));
            
            //  cout << "Merging: " << endl << get<2>(distRepoAB(slpB->root(), slpA->root())) << endl << " with: " << endl
            //       << get<2>(distRepoAC(slpC->root(), slpA_->root())) << endl;
            
            auto merged = Merger::verticalMerge(&distRepoAB(slpB->root(), slpA->root()),
                                                &distRepoAC(slpC->root(), slpA_->root()));
            
            /*
            cout << "Should be: " << endl << distRepoABC(slpBC->root(), slpA__->root()).matrix() << endl;
            cout << "Merged: " << endl << merged->matrix() << endl;
             */
            
            //  cout << "Merged: " << endl << get<2>(*merged) << endl;
            
            if (distRepoABC(slpBC->root(), slpA__->root()) != *merged)
              return false;
          }
        }
      }
      
      return true;
    }
    
    template <class BaseDIST, class TestDIST, class SLPBuilder = SLP::SimpleCompressionSLPBuilder>
    static bool test_dist_repository()
    {
      for (int64_t length = 10; length < 100; length = max((int64_t)((double)length * 1.8), length + 1)) {
        string A = Benchmark::generate_string(length - 3 , { 'a', 'b', 'c' });
        string B = Benchmark::generate_string(length + 4, { 'a', 'b', 'c' });
        
        for (int64_t x = 3; x < min(A.length(), B.length()); x = max((int64_t)((double)x * 4), x + 1)) {
          unique_ptr<SLP::SLP> slpA = SLPBuilder::build(A);
          unique_ptr<SLP::SLP> slpB = SLPBuilder::build(B);
          auto partitionA = SLP::Partitioner::partition(*slpA, x);
          auto partitionB = SLP::Partitioner::partition(*slpB, x);
          
          unique_ptr<SLP::SLP> slpA_ = SLPBuilder::build(A);
          unique_ptr<SLP::SLP> slpB_ = SLPBuilder::build(B);
          auto partitionA_ = SLP::Partitioner::partition(*slpA_, x);
          auto partitionB_ = SLP::Partitioner::partition(*slpB_, x);
          
          TestDIST distRepo(*slpA, *slpB);
          distRepo.build(get<1>(partitionA), get<1>(partitionB));
          
          BaseDIST distRepo_(*slpA_, *slpB_);
          distRepo_.build(get<1>(partitionA_), get<1>(partitionB_));
          
          for (int64_t i = 0; i < get<1>(partitionA).size(); ++i) {
            for (int64_t j = 0; j < get<1>(partitionB).size(); ++j) {
              auto a = get<1>(partitionA)[i];
              auto a_ = get<1>(partitionA_)[i];
              auto b = get<1>(partitionB)[j];
              auto b_ = get<1>(partitionB_)[j];
              
              if (distRepo(a, b) != distRepo_(a_, b_)) {
                return false;
              }
            }
          }
        }
      }
      
      return true;
    }
    
    template <class Aligner>
    static bool test_slp_compression_align() {
      for (int64_t length = 10; length < 100; length = max((int64_t)((double)length * 1.8), length + 1)) {
        string A = Benchmark::generate_string(length - 3);
        string B = Benchmark::generate_string(length + 4);
        
        for (double xfactor = 2; xfactor < 3; ++xfactor) {
          Aligner aligner(A, B, xfactor);
          Simple::EditDistance simple(A, B);
          
          const auto found = aligner.edit_distance();
          const auto should_be = simple.calculate();
          if (found != should_be) {
            return false;
          }
        }
      }
      return true;
    }
    
    static bool test_blow_up_method() {
      for (int64_t length = 10; length < 1000; length = max((int64_t)((double)length * 1.2), length + 1)) {
        const string A = Benchmark::generate_string(length - 3);
        const string B = Benchmark::generate_string(length + 4);
        
        // Compute the actual edit distance
        Simple::EditDistance editDist(A, B);
        const int64_t edit_distance = editDist.calculate();
        
        // Blow up
        string A_blown(2 * A.length(), '$');
        for (int64_t i = 0; i < A.length(); ++i) {
          // A_blown[2 * i] = '$';
          A_blown[2 * i + 1] = A[i];
        }
        
        string B_blown(2 * B.length(), '$');
        for (int64_t j = 0;j < B.length(); ++j) {
          // B_blown[2 * h] = '$';
          B_blown[2 * j + 1] = B[j];
        }
        
        // Compute LCS and normalize
        Simple::LCS LCS(A_blown, B_blown);
        const int64_t lcs = LCS.calculate();
        const int64_t edit_distance_ = (int64_t)(A.length() + B.length()) - lcs;
        
        if (edit_distance != edit_distance_)
          return false;
      }
      
      return true;
    }
    
    static bool test_union_find() {
      typedef Utils::UnionFind<int64_t> UF;
      typedef UF::Element Element;
      
      vector<unique_ptr<Element>> elements;
      for (int64_t i = 1; i <= 10; ++i)
        elements.push_back(unique_ptr<Element>(new Element(i)));
      
      auto get = [&elements](int64_t index) -> Element* {
        return elements[index].get();
      };
      
      if (UF::Find(get(0)) != get(0))
        return false;
      
      UF::Union(get(0), get(1));
      UF::Union(get(2), get(3));
      UF::Union(get(0), get(3));
      
      if (UF::Find(get(2)) != UF::Find(get(0)))
        return false;
      
      if (UF::Find(get(9)) == UF::Find(get(0)))
        return false;
        
      return true;
    }
    
    static bool test_interval_union_find() {
      const int64_t N = 10000;
      Utils::IntervalUnionFind iuf(N);
      for (uint64_t i = 0; i < N; ++i) {
        if (iuf.Find(i) != i)
          return false;
      }
      
      for (uint64_t depth = 1; depth < log2(N); ++depth) {
        for (uint64_t i = 0; i < N - pow(2, depth - 1); i += pow(2, depth)) {
          assert(i + pow(2, depth - 1) < N);
          iuf.Union(i, i + pow(2, depth - 1));
        }
        
        for (uint64_t i = 0; i < N - pow(2, depth - 1); i += pow(2, depth)) {
          for (uint64_t j = 0; j < pow(2, depth) && i + j < N; ++j) {
            if (iuf.Find(i + j) != min((uint64_t)(i + pow(2, depth) - 1), (uint64_t)(N - 1))) {
              iuf.print();
              return false;
            }
          }
        }
      }
      
      return true;
    }
    
    static bool test_lz_factorization() {
      for (uint64_t len = 100; len <= 40000; len *= 1.7) {
      // string str = Benchmark::generate_string(len, { 'a' });
        string str = Benchmark::generate_string(len, { 'a', 'b', 'c' });
      
        auto baseFactors = SLP::LZFactorize::naive_lz_factorize(str);
        auto testFactors = SLP::LZFactorize::lz_factorize(str);
        // auto experimentalFactors = SLP::LZFactorize::experimental_fast_scan_lz_factorize(str);
      
        /*
         cout << "Should be:" << endl;
         for (auto factor : baseFactors) {
         cout << str.substr(get<0>(factor), get<1>(factor) - get<0>(factor) + 1) << endl;
         }
         
         cout << endl << "Was:" << endl;
         for (auto factor : testFactors1) {
         cout << str.substr(get<0>(factor), get<1>(factor) - get<0>(factor) + 1) << endl;
         }
        */
        
        if (baseFactors != testFactors) {
          cout << "Test factors failed! :(" << endl;
          return false;
        }
        
        /*
        if (baseFactors != experimentalFactors) {
          cout << "Experimental factors failed!" << endl;
          return false;
        }
         */
      }
      
      return true;
    }
  };
}
