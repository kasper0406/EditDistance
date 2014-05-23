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
#include "../benchmark/benchmarker.hpp"


namespace Compression {
  using namespace std;
  
  template <class Factory, class SLPCompressor, class DISTRepo>
  class Aligner {
  public:
    Aligner(string A, string B, double xfactor, double ffactor) : xfactor_(xfactor), ffactor_(ffactor), stats(nullptr), A_(A), B_(B)
    { }
    
    virtual ~Aligner() { };
    
    static vector<pair<string, function<int64_t()>>> run(tuple<string, string, double, double> input, Benchmark::Stats* stats = nullptr) {
      Aligner* aligner = Factory::getInstance(input);
      aligner->stats = stats;
      
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
    double ffactor_;
    
    Benchmark::Stats* stats;
    
    string A_, B_;
    unique_ptr<SLP::SLP> slpA_, slpB_;
    vector<SLP::Production*> partitionA_, partitionB_;
    vector<SLP::Production*> blocksA_, blocksB_;
    
    unique_ptr<DISTRepo> distRepo;
  };
  
  template <class SLPCompressor, class DISTRepo>
  class EditDistanceAligner : public Aligner<EditDistanceAligner<SLPCompressor, DISTRepo>, SLPCompressor, DISTRepo> {
  public:
    EditDistanceAligner(string A, string B, double xfactor, double ffactor) : Aligner<EditDistanceAligner, SLPCompressor, DISTRepo>(A, B, xfactor, ffactor)
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
        const int64_t b_len = b->associatedStringLen;
        
        vector<int64_t> block_row(b_len + 1, 0);
        for (int64_t i = 0; i <= b_len; i++)
          block_row[i] = Bpos + i;
        
        int64_t Apos = 0;
        for (auto a : this->partitionA_) {
          const int64_t a_len = a->associatedStringLen;
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
    
    static string short_name() {
      return "editdistalign";
    }
    
    static Aligner<EditDistanceAligner, SLPCompressor, DISTRepo>* getInstance(tuple<string, string, double, double> input) {
      auto instance = new EditDistanceAligner<SLPCompressor, DISTRepo>();
      tie(instance->A_, instance->B_, instance->xfactor_, instance->ffactor_) = input;
      return instance;
    }
    
  private:
    EditDistanceAligner() { }
    
    void generateSLPs() {
      this->slpA_ = move(SLPCompressor::build(this->A_));
      this->slpB_ = move(SLPCompressor::build(this->B_));
      
      cout << "Compression factor A: " << this->slpA_->compressionFactor() << endl;
      cout << "Compression factor B: " << this->slpB_->compressionFactor() << endl;
      
      auto findX = [this] (SLP::SLP* slp) {
        const double f =  max((double)2, this->ffactor_ * ((double)slp->derivedLength() / (double)slp->productions()));
        const double x = max(f / sqrt(log2(f)), min(5., (double)slp->derivedLength()));
        
        return min((int64_t)(x * this->xfactor_), slp->derivedLength());
      };
      
      if (this->stats != nullptr) {
        this->stats->A_derivedLength = this->slpA_->derivedLength();
        this->stats->B_derivedLength = this->slpB_->derivedLength();
        this->stats->A_productions = this->slpA_->productions();
        this->stats->B_productions = this->slpB_->productions();
      }
      
      tie(this->partitionA_, this->blocksA_) = SLP::Partitioner::partition(*this->slpA_, findX(this->slpA_.get()));
      tie(this->partitionB_, this->blocksB_) = SLP::Partitioner::partition(*this->slpB_, findX(this->slpB_.get()));
    }
  };
  
  template <class SLPCompressor, class DISTRepo>
  class LCSBlowUpAligner : public Aligner<LCSBlowUpAligner<SLPCompressor, DISTRepo>, SLPCompressor, DISTRepo> {
  public:
    LCSBlowUpAligner(string A, string B, double xfactor, double ffactor) : Aligner<LCSBlowUpAligner, SLPCompressor, DISTRepo>(A, B, xfactor, ffactor)
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
        const int64_t b_len = b->associatedStringLen;
        vector<int64_t> block_row(b_len + 1, 0);
        
        int64_t Apos = 0;
        for (auto a : this->partitionA_) {
          const int64_t a_len = a->associatedStringLen;
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
          const uint64_t res = (this->slpA_->derivedLength() + this->slpB_->derivedLength()) / 2 - result;
          if (this->stats != nullptr)
            this->stats->result = res;
          return res;
        }
      }
      
      throw runtime_error("This should not happen!");
    }
    
    static string name() {
      return "LCSBlowUpAligner(" + SLPCompressor::name() + ", " + DISTRepo::name() + ")";
    }
    
    static string short_name() {
      return "lcs_blowup";
    }
    
    static Aligner<LCSBlowUpAligner, SLPCompressor, DISTRepo>* getInstance(tuple<string, string, double, double> input) {
      auto instance = new LCSBlowUpAligner<SLPCompressor, DISTRepo>();
      tie(instance->A_, instance->B_, instance->xfactor_, instance->ffactor_) = input;
      return instance;
    }
    
  private:
    LCSBlowUpAligner() { }
    
    void generateSLPs() {
      this->slpA_ = move(SLPCompressor::build(this->A_));
      // cout << "Compression factor A before blow-up: " << this->slpA_->compressionFactor() << endl;
      SLP::BlowUpSLPTransformer::blowUpSLP(this->slpA_.get());
      
      this->slpB_ = move(SLPCompressor::build(this->B_));
      // cout << "Compression factor B before blow-up: " << this->slpB_->compressionFactor() << endl;
      SLP::BlowUpSLPTransformer::blowUpSLP(this->slpB_.get());
      
      /*
      cout << "Compression factor A: " << this->slpA_->compressionFactor() << endl;
      cout << "Compression factor B: " << this->slpB_->compressionFactor() << endl;
       */
      
      auto findX = [this] (SLP::SLP* slp) -> int64_t {
        const double f =  max((double)2, this->ffactor_ * ((double)slp->derivedLength() / (double)slp->productions()));
        const double x = max(f / sqrt(log2(f)), min(5., (double)slp->derivedLength()));
        
        return max((int64_t)1, min((int64_t)(x * this->xfactor_), slp->derivedLength()));
        
        /*
        // return min((int64_t)5, slp->derivedLength());
        
        const int64_t f = max((int64_t)2, slp->derivedLength() / slp->productions());
        const int64_t x = max(f / (int64_t)sqrt(log2(f)), min((int64_t)5, slp->derivedLength()));
        
        return max(int64_t(1), min((int64_t)(x * this->xfactor_), slp->derivedLength()));
         */
      };
      
      const uint64_t A_x = findX(this->slpA_.get());
      const uint64_t A_y = findX(this->slpB_.get());
      
      if (this->stats != nullptr) {
        this->stats->A_derivedLength = this->slpA_->derivedLength();
        this->stats->B_derivedLength = this->slpB_->derivedLength();
        this->stats->A_productions = this->slpA_->productions();
        this->stats->B_productions = this->slpB_->productions();
        this->stats->A_x = A_x;
        this->stats->A_y = A_y;
      }
      
      /*
      cout << "A x: " << findX(this->slpA_.get()) << endl;
      cout << "B x: " << findX(this->slpB_.get()) << endl;
       */
      
      tie(this->partitionA_, this->blocksA_) = SLP::Partitioner::partition(*this->slpA_, A_x);
      tie(this->partitionB_, this->blocksB_) = SLP::Partitioner::partition(*this->slpB_, A_y);
    }
  };
}
