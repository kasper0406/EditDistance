#pragma once

#include <vector>
#include <cstdint>
#include <tuple>
#include <string>
#include <stdexcept>
#include <memory>
#include <utility>
#include <list>

#include "SLP.hpp"
#include "DIST.hpp"

namespace Compression {
  using namespace std;
  
  template <class Factory, class SLPCompressor, class DISTRepo>
  class Aligner {
  public:
    Aligner(string A, string B, int64_t x) : x_(x), A_(A), B_(B)
    {
      generateSLPs();
      buildDISTRepo();
    }
    
    virtual ~Aligner() { };
    
    static vector<pair<string, function<int64_t()>>> run(tuple<string, string, int64_t> input) {
      Aligner* aligner = Factory::getInstance(input);
      
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
    
    virtual int64_t edit_distance() = 0;
    
  protected:
    Aligner() { }
    
    virtual void generateSLPs() {
      slpA_ = move(SLPCompressor::build(A_));
      tie(partitionA_, blocksA_) = SLP::Partitioner::partition(*slpA_, x_);
      
      slpB_ = move(SLPCompressor::build(B_));
      tie(partitionB_, blocksB_) = SLP::Partitioner::partition(*slpB_, x_);
    }
    
    void buildDISTRepo() {
      distRepo = unique_ptr<DISTRepo>(new DISTRepo(*slpA_, *slpB_));
      distRepo->build(blocksA_, blocksB_);
    }
    
    int64_t x_;
    
    string A_, B_;
    unique_ptr<SLP::SLP> slpA_, slpB_;
    vector<SLP::Production*> partitionA_, partitionB_;
    vector<SLP::Production*> blocksA_, blocksB_;
    
    unique_ptr<DISTRepo> distRepo;
  };
  
  template <class SLPCompressor, class DISTRepo>
  class EditDistanceAligner : public Aligner<EditDistanceAligner<SLPCompressor, DISTRepo>, SLPCompressor, DISTRepo> {
  public:
    EditDistanceAligner(string A, string B, int64_t x) : Aligner<EditDistanceAligner, SLPCompressor, DISTRepo>(A, B, x) { }
    
    int64_t edit_distance() {
      // Base case for column
      vector<int64_t> prev_column(this->A_.size() + 1, 0), cur_column(this->A_.size() + 1, 0);
      for (int64_t Apos = 0; Apos <= this->A_.size(); ++Apos) {
        prev_column[Apos] = Apos;
      }
      
      int64_t Bpos = 0;
      for (auto b : this->partitionB_) {
        const int64_t b_len = b->associatedString.size();
        
        vector<int64_t> block_row(b_len + 1, 0);
        for (int64_t i = 0; i <= b_len; i++)
          block_row[i] = Bpos + i;
        
        int64_t Apos = 0;
        for (auto a : this->partitionA_) {
          const int64_t a_len = a->associatedString.size();
          assert(block_row[0] == prev_column[Apos]);
          
          // Build inputs to block
          vector<int64_t> I(a_len + b_len + 1, 0);
          for (int64_t i = 0; i <= a_len; ++i) {
            I[i] = prev_column[Apos + a_len - i];
            // cout << "I[" << i << "] = " << I[i] << endl;
          }
          for (int64_t i = a_len + 1; i < I.size(); ++i)
            I[i] = block_row[i - a_len];
          
          // Computes outputs from block
          auto O = this->distRepo->apply(a, b, I);
          
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
        if (Bpos == this->B_.size())
          return block_row[block_row.size() - 1];
      }
      
      throw runtime_error("This should not happen!");
    }
    
    static string name() {
      return "SLPAlign(" + SLPCompressor::name() + ", " + DISTRepo::name() + ")";
    }
    
    static Aligner<EditDistanceAligner, SLPCompressor, DISTRepo>* getInstance(tuple<string, string, int64_t> input) {
      auto instance = new EditDistanceAligner<SLPCompressor, DISTRepo>();
      tie(instance->A_, instance->B_, instance->x_) = input;
      return instance;
    }
    
  private:
    EditDistanceAligner() { }
  };
  
  template <class SLPCompressor, class DISTRepo>
  class LCSBlowUpAligner : public Aligner<LCSBlowUpAligner<SLPCompressor, DISTRepo>, SLPCompressor, DISTRepo> {
  public:
    /**
     * Consider making the blow up directly on the SLP potentially saving some time.
     */
    LCSBlowUpAligner(string A, string B, int64_t x) : Aligner<LCSBlowUpAligner, SLPCompressor, DISTRepo>(blowUp(A), blowUp(B), x) { }
    
    int64_t edit_distance() {
      // Base case for column
      vector<int64_t> prev_column(this->A_.size() + 1, 0), cur_column(this->A_.size() + 1, 0);
      
      int64_t Bpos = 0;
      for (auto b : this->partitionB_) {
        const int64_t b_len = b->associatedString.size();
        vector<int64_t> block_row(b_len + 1, 0);
        
        int64_t Apos = 0;
        for (auto a : this->partitionA_) {
          const int64_t a_len = a->associatedString.size();
          assert(block_row[0] == prev_column[Apos]);
          
          // Build inputs to block
          vector<int64_t> I(a_len + b_len + 1, 0);
          for (int64_t i = 0; i <= a_len; ++i) {
            I[i] = prev_column[Apos + a_len - i];
            // cout << "I[" << i << "] = " << I[i] << endl;
          }
          for (int64_t i = a_len + 1; i < I.size(); ++i)
            I[i] = block_row[i - a_len];
          
          // Computes outputs from block
          auto O = this->distRepo->apply(a, b, I);
          
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
        if (Bpos == this->B_.size()) {
          const int64_t result = block_row[block_row.size() - 1];
          return (this->A_.length() + this->B_.length()) / 2 - result;
        }
      }
      
      throw runtime_error("This should not happen!");
    }
    
    static string name() {
      return "LCSBlowUpAligner(" + SLPCompressor::name() + ", " + DISTRepo::name() + ")";
    }
    
    static string blowUp(string S) {
      string result(2 * S.length(), '$');
      for (int64_t i = 0; i < S.length(); ++i) {
        // result[2 * i] = '$';
        result[2 * i + 1] = S[i];
      }
      return result;
    }
    
    static Aligner<LCSBlowUpAligner, SLPCompressor, DISTRepo>* getInstance(tuple<string, string, int64_t> input) {
      auto instance = new LCSBlowUpAligner<SLPCompressor, DISTRepo>();
      instance->A_ = blowUp(get<0>(input));
      instance->B_ = blowUp(get<1>(input));
      instance->x_ = get<2>(input);
      return instance;
    }
    
  private:
    LCSBlowUpAligner() { }
  };
}
