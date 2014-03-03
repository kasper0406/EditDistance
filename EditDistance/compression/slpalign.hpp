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
    Aligner(string A, string B, double xfactor) : xfactor_(xfactor), A_(A), B_(B)
    { }
    
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
    
    virtual void generateSLPs() = 0;
    
    void buildDISTRepo() {
      distRepo = unique_ptr<DISTRepo>(new DISTRepo(*slpA_, *slpB_));
      distRepo->build(blocksA_, blocksB_);
    }
    
    double xfactor_;
    
    string A_, B_;
    unique_ptr<SLP::SLP> slpA_, slpB_;
    vector<SLP::Production*> partitionA_, partitionB_;
    vector<SLP::Production*> blocksA_, blocksB_;
    
    unique_ptr<DISTRepo> distRepo;
  };
  
  template <class SLPCompressor, class DISTRepo>
  class EditDistanceAligner : public Aligner<EditDistanceAligner<SLPCompressor, DISTRepo>, SLPCompressor, DISTRepo> {
  public:
    EditDistanceAligner(string A, string B, double xfactor) : Aligner<EditDistanceAligner, SLPCompressor, DISTRepo>(A, B, xfactor)
    {
      this->generateSLPs();
      this->buildDISTRepo();
    }
    
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
    
    static Aligner<EditDistanceAligner, SLPCompressor, DISTRepo>* getInstance(tuple<string, string, double> input) {
      auto instance = new EditDistanceAligner<SLPCompressor, DISTRepo>();
      tie(instance->A_, instance->B_, instance->xfactor_) = input;
      return instance;
    }
    
  private:
    EditDistanceAligner() { }
    
    void generateSLPs() {
      this->slpA_ = move(SLPCompressor::build(this->A_));
      this->slpB_ = move(SLPCompressor::build(this->B_));
      
      auto findX = [this] (SLP::SLP* slp) {
        const int64_t f = slp->derivedLength() / slp->productions();
        const int64_t x = max(f / (int64_t)sqrt(log2(f)), min((int64_t)5, slp->derivedLength()));
        
        cout << "x: " << (int64_t)(x * this->xfactor_) << endl;
        
        return min((int64_t)(x * this->xfactor_), slp->derivedLength());
      };
      
      tie(this->partitionA_, this->blocksA_) = SLP::Partitioner::partition(*this->slpA_, findX(this->slpA_.get()));
      tie(this->partitionB_, this->blocksB_) = SLP::Partitioner::partition(*this->slpB_, findX(this->slpB_.get()));
    }
  };
  
  template <class SLPCompressor, class DISTRepo>
  class LCSBlowUpAligner : public Aligner<LCSBlowUpAligner<SLPCompressor, DISTRepo>, SLPCompressor, DISTRepo> {
  public:
    LCSBlowUpAligner(string A, string B, double xfactor) : Aligner<LCSBlowUpAligner, SLPCompressor, DISTRepo>(A, B, xfactor)
    {
      this->generateSLPs();
      this->buildDISTRepo();
    }
    
    int64_t edit_distance() {
      assert(2 * this->A_.size() == this->slpA_->derivedLength() && 2 * this->B_.size() == this->slpB_->derivedLength());
      
      // Base case for column
      vector<int64_t> prev_column(this->slpA_->derivedLength() + 1, 0), cur_column(this->slpA_->derivedLength() + 1, 0);
      
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
        if (Bpos == this->slpB_->derivedLength()) {
          const int64_t result = block_row[block_row.size() - 1];
          return (this->slpA_->derivedLength() + this->slpB_->derivedLength()) / 2 - result;
        }
      }
      
      throw runtime_error("This should not happen!");
    }
    
    static string name() {
      return "LCSBlowUpAligner(" + SLPCompressor::name() + ", " + DISTRepo::name() + ")";
    }
    
    static Aligner<LCSBlowUpAligner, SLPCompressor, DISTRepo>* getInstance(tuple<string, string, double> input) {
      auto instance = new LCSBlowUpAligner<SLPCompressor, DISTRepo>();
      instance->A_ = get<0>(input);
      instance->B_ = get<1>(input);
      instance->xfactor_ = get<2>(input);
      return instance;
    }
    
  private:
    LCSBlowUpAligner() { }
    
    void generateSLPs() {
      this->slpA_ = move(SLPCompressor::build(this->A_));
      SLP::BlowUpSLPTransformer::blowUpSLP(this->slpA_.get());
      
      this->slpB_ = move(SLPCompressor::build(this->B_));
      SLP::BlowUpSLPTransformer::blowUpSLP(this->slpB_.get());
      
      auto findX = [this] (SLP::SLP* slp) {
        const int64_t f = slp->derivedLength() / slp->productions();
        const int64_t x = max(f / (int64_t)sqrt(log2(f)), min((int64_t)5, slp->derivedLength()));
        
        cout << "x: " << (int64_t)(x * this->xfactor_) << endl;
        
        return min((int64_t)(x * this->xfactor_), slp->derivedLength());
      };
      
      tie(this->partitionA_, this->blocksA_) = SLP::Partitioner::partition(*this->slpA_, findX(this->slpA_.get()));
      tie(this->partitionB_, this->blocksB_) = SLP::Partitioner::partition(*this->slpB_, findX(this->slpB_.get()));
    }
  };
}
