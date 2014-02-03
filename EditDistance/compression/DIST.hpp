#pragma once

#include <vector>
#include <tuple>
#include <list>
#include <memory>

#include "../utils/matrix.hpp"
#include "SLP.hpp"

namespace Compression {
  namespace DIST {
    typedef SLP::SLP StraightLineProgram;
    
    typedef tuple<int64_t, int64_t, Matrix<uint64_t>> DISTTable;
    
    using namespace std;
    using namespace SLP;
    
    class SimpleDISTRepository {
    public:
      SimpleDISTRepository(const StraightLineProgram& slpA,
                           const StraightLineProgram& slpB)
      { }
      
      /**
       * Constructs DIST repository given blocks in both strings.
       */
      void build(const vector<Production*>& A,
                 const vector<Production*>& B) {
        repo_.reserve(A.size());
        
        for (auto a : A) {
          vector<DISTTable> row;
          row.reserve(B.size());
          assert(a->DISTTableIndex == -1);
          a->DISTTableIndex = repo_.size();
          
          for (auto& b : B) {
            assert(b->DISTTableIndex == -1 || row.size() == b->DISTTableIndex);
            b->DISTTableIndex = row.size();
            
            const int64_t m = a->associatedString.size();
            const int64_t n = b->associatedString.size();
            Matrix<uint64_t> matrix(m + n + 1, m + n + 1, numeric_limits<uint64_t>::max());
            
            // Fill up the DIST matrix
            for (int64_t in = 0; in <= m + n; in++) {
              for (int64_t out = 0; out <= m + n; out++) {
                const int64_t a_start = max(m - in, (int64_t)0);
                const int64_t a_stop = min((m + n) - out, m);
                if (a_start > a_stop) continue; // Invalid substring
                
                const int64_t b_start = max(in - (m), (int64_t)0);
                const int64_t b_stop = min(out, n);
                if (b_start > b_stop) continue; // Invalid substring
                
                /*
                 cout << "(in, out) = (" << in << ", " << out << ")" << endl;
                 cout << a_start << " -> " << a_stop << ": " << a->associatedString.substr(a_start, a_stop - a_start) << endl;
                 cout << b_start << " -> " << b_stop << ": " << b->associatedString.substr(b_start, b_stop - b_start) << endl;
                 */
                
                Simple calculator(a->associatedString.substr(a_start, a_stop - a_start),
                                  b->associatedString.substr(b_start, b_stop - b_start));
                // cout << "Res: " << calculator.edit_distance() << endl;
                matrix(in, out) = calculator.edit_distance();
              }
            }
            
            /*
             cout << a->associatedString << endl << b->associatedString << endl;
             cout << matrix << endl << endl;
             */
            
            DISTTable table = { m, n, move(matrix) };
            row.push_back(move(table));
          }
          
          repo_.push_back(move(row));
        }
      }
      
      /**
       * Given two productions A, B and input I, return the output of applying the DIST table identified to A, B to I.
       */
      vector<uint64_t> apply(Production* A, Production* B, const vector<uint64_t>& I) {
        assert(A->DISTTableIndex != -1 && B->DISTTableIndex != -1);
        assert(A->associatedString.size() + B->associatedString.size() + 1 == I.size());
        
        /*
         cout << "A: " << A->associatedString << endl << "B: " << B->associatedString << endl;
         cout << "Associated DIST:" << endl << repo_[A->DISTTableIndex][B->DISTTableIndex] << endl << endl;
         */
        
        vector<uint64_t> O;
        O.reserve(I.size());
        
        for (uint64_t j = 0; j < I.size(); j++) {
          uint64_t minimum = numeric_limits<uint64_t>::max();
          for (uint64_t i = 0; i < I.size(); i++) {
            /*
             cout << "out: " << j << endl;
             cout << "in: " << i << endl;
             cout << "Input: " << I[i] << endl;
             cout << "DIST: " << repo_[A->DISTTableIndex][B->DISTTableIndex](i, j) << endl;
             */
            
            if (get<2>(repo_[A->DISTTableIndex][B->DISTTableIndex])(i, j) == numeric_limits<uint64_t>::max()) continue;
            minimum = min(minimum, I[i] + get<2>(repo_[A->DISTTableIndex][B->DISTTableIndex])(i, j));
          }
          O.push_back(minimum);
        }
        
        return O;
      }
      
      const DISTTable& operator()(Production* a, Production* b) const {
        const int64_t i = a->DISTTableIndex;
        const int64_t j = b->DISTTableIndex;
        
        assert(i != -1 && j != -1 && i < repo_.size() && j < repo_[i].size());
        return repo_[i][j];
      }
      
    private:
      vector<vector<DISTTable>> repo_;
    };
    
    template <class Merger>
    class MergingDISTRepository {
    public:
      MergingDISTRepository(const StraightLineProgram& slpA,
                            const StraightLineProgram& slpB)
        : a_dist_index_count_(0), b_dist_index_count_(0),
          repo_(slpA.productions(), slpB.productions(), nullptr)
      { }
      
      MergingDISTRepository(const MergingDISTRepository&& other)
        : repo_(move(other.repo_))
      { }
      
      MergingDISTRepository(const MergingDISTRepository& other)
        : repo_(0, 0, nullptr)
      {
        throw runtime_error("Do not copy a DIST repo!");
      }
      
      MergingDISTRepository& operator=(const MergingDISTRepository&& other) {
        repo_ = move(other.repo_);
        return *this;
      }
      
      MergingDISTRepository& operator=(const MergingDISTRepository& other) {
        throw runtime_error("Do not copy a DIST repo!");
      }
      
      void build(const vector<Production*>& A,
                 const vector<Production*>& B) {
        // Handle case where productions derive exactly the associated string
        for (auto a : A) {
          for (auto b : B) {
            build(a, b);
          }
        }
      }
      
      /**
       * Given two productions A, B and input I, return the output of applying the DIST table identified to A, B to I.
       */
      vector<uint64_t> apply(Production* A, Production* B, const vector<uint64_t>& I) {
        // TODO: Factor this code out into the DIST table itself.
        assert(A->DISTTableIndex != -1 && B->DISTTableIndex != -1);
        assert(A->associatedString.size() + B->associatedString.size() + 1 == I.size());
        
        /*
         cout << "A: " << A->associatedString << endl << "B: " << B->associatedString << endl;
         cout << "Associated DIST:" << endl << repo_[A->DISTTableIndex][B->DISTTableIndex] << endl << endl;
         */
        
        vector<uint64_t> O;
        O.reserve(I.size());
        
        for (uint64_t j = 0; j < I.size(); j++) {
          uint64_t minimum = numeric_limits<uint64_t>::max();
          for (uint64_t i = 0; i < I.size(); i++) {
            /*
             cout << "out: " << j << endl;
             cout << "in: " << i << endl;
             cout << "Input: " << I[i] << endl;
             cout << "DIST: " << repo_[A->DISTTableIndex][B->DISTTableIndex](i, j) << endl;
             */
            
            assert(repo_(A->DISTTableIndex, B->DISTTableIndex) != nullptr);
            if (get<2>(*repo_(A->DISTTableIndex, B->DISTTableIndex))(i, j) == numeric_limits<uint64_t>::max()) continue;
            minimum = min(minimum, I[i] + get<2>(*repo_(A->DISTTableIndex, B->DISTTableIndex))(i, j));
          }
          O.push_back(minimum);
        }
        
        return O;
      }
      
      const DISTTable& operator()(Production* a, Production* b) const {
        const int64_t i = a->DISTTableIndex;
        const int64_t j = b->DISTTableIndex;
        
        assert(i != -1 && j != -1 && i < repo_.getRows() && j < repo_.getColumns() && repo_(i, j) != nullptr);
        return *repo_(i, j);
      }
      
    private:
      shared_ptr<DISTTable> build(Terminal* a, Terminal* b)
      {
        // cout << "Term - Term" << endl;
        
        const uint64_t inf = numeric_limits<uint64_t>::max();
        const uint64_t is_equal = a->symbol() == b->symbol();
        auto result = shared_ptr<DISTTable>(new DISTTable(1, 1, Matrix<uint64_t>(3, 3, {
          { 0, 1, inf },
          { 1, (1 - is_equal), 1 },
          { inf, 1, 0 }
        })));
        repo_(a->DISTTableIndex, b->DISTTableIndex) = result;
        return result;
      }
      
      shared_ptr<DISTTable> build(NonTerminal* a, const Terminal* b)
      {
        // cout << "NonTerm - Term" << endl;
        const bool a_type1 = a->associatedString.size() == a->derivedStringLength;
        
        Production *a_left, *a_right;
        if (a_type1) {
          a_left = a->left();
          a_right = a->right();
        } else {
          // TODO: Consider if this can happen!
          cout << "COVERED 1!" << endl;
          assert(!a->type != KEY);
          if (a->type == LEFT) {
            if (a->TYPE2_next == nullptr)
              return repo_(a->DISTTableIndex, b->DISTTableIndex) = build(a->right(), (Production*)b);
            
            a_left = a->TYPE2_next;
            a_right = a->right();
          } else { // RIGHT
            if (a->TYPE2_next == nullptr)
              return repo_(a->DISTTableIndex, b->DISTTableIndex) = build(a->left(), (Production*)b);
            
            a_left = a->left();
            a_right = a->TYPE2_next;
          }
        }
        
        auto top = build(a_left, (Production*)b);
        auto bottom = build(a_right, (Production*)b);
        auto dist = Merger::verticalMerge(top.get(), bottom.get());
        return repo_(a->DISTTableIndex, b->DISTTableIndex) = dist;
      }
      
      shared_ptr<DISTTable> build(Terminal* a, NonTerminal* b)
      {
        // cout << "Term - NonTerm" << endl;
        const bool b_type1 = b->associatedString.size() == b->derivedStringLength;
        
        Production *b_left, *b_right;
        if (b_type1) {
          b_left = b->left();
          b_right = b->right();
        } else {
          // TODO: Consider if this can happen!
          cout << "COVERED 2!" << endl;
          assert(b->type != KEY);
          if (b->type == LEFT) {
            if (b->TYPE2_next == nullptr)
              return repo_(a->DISTTableIndex, b->DISTTableIndex) = build(a, b->right());
            
            b_left = b->TYPE2_next;
            b_right = b->right();
          } else { // RIGHT
            if (b->TYPE2_next == nullptr)
              return repo_(a->DISTTableIndex, b->DISTTableIndex) = build(a, b->left());
            
            b_left = b->left();
            b_right = b->TYPE2_next;
          }
        }
        
        auto left = build((Production*)a, b_left);
        auto right = build((Production*)a, b_right);
        auto dist = Merger::horizontalMerge(left.get(), right.get());
        return repo_(a->DISTTableIndex, b->DISTTableIndex) = dist;
      }
      
      shared_ptr<DISTTable> build(NonTerminal* a, NonTerminal* b)
      {
        // cout << "NonTerm - NonTerm" << endl;
        
        const bool a_type1 = a->associatedString.size() == a->derivedStringLength;
        const bool b_type1 = b->associatedString.size() == b->derivedStringLength;
        
        // Determine things to merge
        Production *a_left, *a_right, *b_left, *b_right;
        if (a_type1) {
          a_left = a->left();
          a_right = a->right();
        } else {
          assert(a->type != KEY);
          if (a->type == LEFT) {
            if (a->TYPE2_next == nullptr)
              return repo_(a->DISTTableIndex, b->DISTTableIndex) = build(a->right(), b);
            
            a_left = a->TYPE2_next;
            a_right = a->right();
          } else { // RIGHT
            if (a->TYPE2_next == nullptr)
              return repo_(a->DISTTableIndex, b->DISTTableIndex) = build(a->left(), b);
            
            a_left = a->left();
            a_right = a->TYPE2_next;
          }
        }
        
        if (b_type1) {
          b_left = b->left();
          b_right = b->right();
        } else {
          assert(b->type != KEY);
          if (b->type == LEFT) {
            if (b->TYPE2_next == nullptr)
              return repo_(a->DISTTableIndex, b->DISTTableIndex) = build(a, b->right());
            
            b_left = b->TYPE2_next;
            b_right = b->right();
          } else { // RIGHT
            if (b->TYPE2_next == nullptr)
              return repo_(a->DISTTableIndex, b->DISTTableIndex) = build(a, b->left());
            
            b_left = b->left();
            b_right = b->TYPE2_next;
          }
        }
        
        // Recursively build DIST combinations of children
        auto ll = build(a_left, b_left);
        auto lr = build(a_left, b_right);
        auto rl = build(a_right, b_left);
        auto rr = build(a_right, b_right);
        
        // Merge
        auto top = Merger::horizontalMerge(ll.get(), lr.get());
        auto bottom = Merger::horizontalMerge(rl.get(), rr.get());
        auto dist = Merger::verticalMerge(top.get(), bottom.get());
          
        return repo_(a->DISTTableIndex, b->DISTTableIndex) = dist;
      }
      
      shared_ptr<DISTTable> build(Production* a, Production* b)
      {
        // cout << "'Double dispatch': ";
        
        if (a->DISTTableIndex == -1) a->DISTTableIndex = a_dist_index_count_++;
        if (b->DISTTableIndex == -1) b->DISTTableIndex = b_dist_index_count_++;
        if (repo_(a->DISTTableIndex, b->DISTTableIndex) != nullptr) { // Already computed
          // cout << "No dispatch - already computed :)." << endl;
          return repo_(a->DISTTableIndex, b->DISTTableIndex);
        }
        
        // Dispatch to correct method.
        // TODO: Consider if this can be simplified! This is pretty ugly right now.
        NonTerminal* a_nonterm = dynamic_cast<NonTerminal*>(a);
        NonTerminal* b_nonterm = dynamic_cast<NonTerminal*>(b);
        
        if (a_nonterm == nullptr) {
          if (b_nonterm == nullptr) {
            return build(static_cast<Terminal*>(a), static_cast<Terminal*>(b));
          } else {
            return build(static_cast<Terminal*>(a), b_nonterm);
          }
        } else {
          if (b_nonterm == nullptr) {
            return build(a_nonterm, static_cast<Terminal*>(b));
          } else {
            return build(a_nonterm, b_nonterm);
          }
        }
      }
      
      uint64_t a_dist_index_count_, b_dist_index_count_;
      Matrix<shared_ptr<DISTTable>> repo_;
    };
    
    class SimpleMerger {
    public:
      static shared_ptr<DISTTable> horizontalMerge(const DISTTable* left, const DISTTable* right) {
        assert(get<0>(*left) == get<0>(*right));
        const int64_t rows = get<0>(*left);
        const int64_t cols = get<1>(*left) + get<1>(*right);
        const int64_t N = rows + cols + 1;
        
        auto result = shared_ptr<DISTTable>(new DISTTable({ rows, cols, Matrix<uint64_t>(N, N, 9) }));
        for (int64_t in = 0; in < N; ++in) {
          for (int64_t out = 0; out < N; ++out) {
            // Only in left table
            if (out <= get<1>(*left)) {
              if (in < get<2>(*left).getRows())
                get<2>(*result)(in, out) = get<2>(*left)(in, out);
              else
                get<2>(*result)(in, out) = numeric_limits<uint64_t>::max();
            } else if (in < get<2>(*left).getRows()) {
              // Spanning both tables
              int64_t left_in = in;
              int64_t right_out = out - get<1>(*left);
              
              uint64_t minimum = numeric_limits<uint64_t>::max();
              for (int64_t right_in = get<0>(*right); right_in >= max(right_out - get<1>(*right), (int64_t)0); --right_in) {
                uint64_t left_out = right_in + get<1>(*left);
                
                uint64_t left_path = get<2>(*left)(left_in, left_out);
                uint64_t right_path = get<2>(*right)(right_in, right_out);
                if (left_path == numeric_limits<uint64_t>::max() ||
                    right_path == numeric_limits<uint64_t>::max())
                {
                  continue;
                }
                
                minimum = min(minimum, left_path + right_path);
              }
              
              get<2>(*result)(in, out) = minimum;
            } else {
              // Only in right table
              const int64_t right_in = in - get<2>(*left).getRows() + get<0>(*right) + 1;
              const int64_t right_out = out - get<1>(*left);
              get<2>(*result)(in, out) = get<2>(*right)(right_in, right_out);
            }
          }
        }
        
        return result;
      }
      
      static shared_ptr<DISTTable> verticalMerge(const DISTTable* top, const DISTTable* bottom) {
        assert(get<1>(*top) == get<1>(*bottom));
        const int64_t rows = get<0>(*bottom) + get<0>(*top);
        const int64_t cols = get<1>(*bottom);
        const int64_t N = rows + cols + 1;
        
        auto result = shared_ptr<DISTTable>(new DISTTable({ rows, cols, Matrix<uint64_t>(N, N, 9) }));
        for (int64_t in = 0; in < N; ++in) {
          for (int64_t out = 0; out < N; ++out) {
            if (out >= get<2>(*bottom).getRows()) {
              // Only in top table
              if (in > get<0>(*bottom)) {
                const int64_t top_in = in - get<0>(*bottom);
                const int64_t top_out = out - get<0>(*bottom);
                get<2>(*result)(in, out) = get<2>(*top)(top_in, top_out);
              } else {
                get<2>(*result)(in, out) = numeric_limits<uint64_t>::max();
              }
            } else if (in > get<0>(*bottom)) {
              // Spanning both
              uint64_t minimum = numeric_limits<uint64_t>::max();
              for (int64_t k = 0; k <= min(out, get<1>(*bottom)); ++k) {
                const uint64_t top_path = get<2>(*top)(in - get<0>(*bottom), k);
                const uint64_t bottom_path = get<2>(*bottom)(get<0>(*bottom) + k, out);
                if (top_path == numeric_limits<uint64_t>::max() ||
                    bottom_path == numeric_limits<uint64_t>::max())
                {
                  continue;
                }
                
                minimum = min(minimum, top_path + bottom_path);
              }
              get<2>(*result)(in, out) = minimum;
            } else {
              // Only in bottom table
              get<2>(*result)(in, out) = get<2>(*bottom)(in, out);
            }
          }
        }
        
        return result;
      }
    };
  }
}
