#pragma once

#include <vector>
#include <cstdint>
#include <tuple>
#include <string>
#include <stdexcept>
#include <memory>
#include <utility>

#include "SLP.hpp"
#include "DIST.hpp"

namespace Compression {
  using namespace std;
  
  template <class SLPCompressor, class DISTRepo>
  class SLPAlign {
  public:
    SLPAlign(string A, string B, uint64_t x) : x_(x), A_(A), B_(B) {
      generateSLPs();
      buildDISTRepo();
    }
    
    static vector<pair<string, function<uint64_t()>>> run(tuple<string, string, uint64_t> input) {
      SLPAlign* aligner = new SLPAlign();
      tie(aligner->A_, aligner->B_, aligner->x_) = input;
      
      auto slp = [aligner]() {
        aligner->generateSLPs();
        return 0;
      };
      
      auto dist = [aligner]() {
        aligner->buildDISTRepo();
        return 0;
      };
      
      auto grid = [aligner]() {
        return aligner->edit_distance();
      };
      
      auto cleanup = [aligner]() {
        delete aligner;
        return 0;
      };
      
      return {
        { "SLP", slp },
        { "DIST", dist },
        { "grid", grid },
        { "cleanup", cleanup }
      };
    }
    
    uint64_t edit_distance() {
      // Base case for column
      vector<uint64_t> prev_column(A_.size() + 1, 0), cur_column(A_.size() + 1, 0);
      for (uint64_t Apos = 0; Apos <= A_.size(); ++Apos) {
        prev_column[Apos] = Apos;
      }
      
      int64_t Bpos = 0;
      for (auto b : partitionB_) {
        const int64_t b_len = b->associatedString.size();
        
        vector<uint64_t> block_row(b_len + 1, 0);
        for (int64_t i = 0; i <= b_len; i++)
          block_row[i] = Bpos + i;
        
        int64_t Apos = 0;
        for (auto a : partitionA_) {
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
          auto O = distRepo.apply(a, b, I);
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
        if (Bpos == B_.size())
          return block_row[block_row.size() - 1];
      }
      
      throw runtime_error("This should not happen!");
    }
    
    static string name() {
      return "SLPAlign";
    }
    
  private:
    SLPAlign() { }
    
    void generateSLPs() {
      slpA_ = move(SLPCompressor::build(A_));
      tie(partitionA_, blocksA_) = SLP::Partitioner::partition(*slpA_, x_);
      
      slpB_ = move(SLPCompressor::build(B_));
      tie(partitionB_, blocksB_) = SLP::Partitioner::partition(*slpB_, x_);
    }
    
    void buildDISTRepo() {
      distRepo.build(blocksA_, blocksB_);
    }
    
    uint64_t x_;
    
    string A_, B_;
    unique_ptr<SLP::SLP> slpA_, slpB_;
    vector<SLP::Production*> partitionA_, blocksA_;
    vector<SLP::Production*> partitionB_, blocksB_;
    
    DISTRepo distRepo;
  };
}
