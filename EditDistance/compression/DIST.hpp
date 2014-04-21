#pragma once

#include <vector>
#include <tuple>
#include <list>
#include <iterator>
#include <memory>

#include "../simple/simple.hpp"

#include "../utils/matrix.hpp"
#include "../utils/math.hpp"
#include "../utils/unionfind.hpp"
#include "SLP.hpp"
#include "maxmultiply.hpp"

namespace Compression {
  namespace DIST {
    typedef ::Compression::SLP::SLP StraightLineProgram;
    using namespace std;
    using namespace ::Compression::SLP;
    using namespace ::Utils::Math;
    
    using namespace ::Simple;
    
    class EditDistanceDISTTable {
    public:
      EditDistanceDISTTable(int64_t rows, int64_t cols, Matrix<int64_t> matrix)
        : rows_(rows), cols_(cols), matrix_(move(matrix))
      { }
      
      int64_t rows() const { return rows_; }
      int64_t cols() const { return cols_; }
      
      int64_t size() const { return rows() + cols(); }
      
      const Matrix<int64_t>& matrix() const { return matrix_; }
      Matrix<int64_t>& matrix() { return matrix_; }
      
      static unique_ptr<EditDistanceDISTTable> BaseCase(bool is_equal) {
        const int64_t inf = numeric_limits<int64_t>::max();
        return unique_ptr<EditDistanceDISTTable>(new EditDistanceDISTTable(1, 1, Matrix<int64_t>(3, 3, {
          { 0, 1, inf },
          { 1, (1 - is_equal), 1 },
          { inf, 1, 0 }
        })));
      }
      
      bool operator==(const EditDistanceDISTTable& other) const {
        return rows_ == other.rows_ && cols_ == other.cols_ && matrix_ == other.matrix_;
      }
      
      bool operator!=(const EditDistanceDISTTable& other) const {
        return !((*this) == other);
      }
      
      vector<int64_t> apply(const vector<int64_t>& I) const {
        vector<int64_t> O;
        O.reserve(I.size());
        
        for (int64_t j = 0; j < I.size(); j++) {
          int64_t minimum = numeric_limits<int64_t>::max();
          for (int64_t i = 0; i < I.size(); i++) {
            if (matrix()(i, j) == numeric_limits<int64_t>::max()) continue;
            minimum = min(minimum, I[i] + matrix()(i, j));
          }
          O.push_back(minimum);
        }
        
        return O;
      }
      
      static string name() {
        return "EditDistanceDIST";
      }
      
    private:
      int64_t rows_;
      int64_t cols_;
      Matrix<int64_t> matrix_;
    };
    
    class SimpleLCSDISTTable {
    public:
      SimpleLCSDISTTable(int64_t rows, int64_t cols, Matrix<int64_t> matrix)
        : rows_(rows), cols_(cols), matrix_(move(matrix))
      { }
      
      int64_t rows() const { return rows_; }
      int64_t cols() const { return cols_; }
      
      int64_t size() const { return rows() + cols(); }
      
      const Matrix<int64_t>& matrix() const { return matrix_; }
      Matrix<int64_t>& matrix() { return matrix_; }
      
      static unique_ptr<SimpleLCSDISTTable> BaseCase(bool is_equal) {
        return unique_ptr<SimpleLCSDISTTable>(new SimpleLCSDISTTable(1, 1, Matrix<int64_t>(3, 3, {
          { 1, 1, 1 },
          { 0, is_equal, 1 },
          { -1, 0, 1 }
        })));
      }
      
      bool operator==(const SimpleLCSDISTTable& other) const {
        return rows_ == other.rows_ && cols_ == other.cols_ && matrix_ == other.matrix_;
      }
      
      bool operator!=(const SimpleLCSDISTTable& other) const {
        return !((*this) == other);
      }
      
      vector<int64_t> apply(const vector<int64_t>& I) const {
        vector<int64_t> O; O.reserve(I.size());
        
        const int64_t x = I.size();
        const int64_t m = rows_;
        const int64_t n = cols_;
        
        // Find maximum path for all different paths between inputs and outputs.
        {
          const auto& H = matrix();
          Matrix<int64_t> input(1, x, { I });
          auto bottom_left = Multiplier::multiply(input, H,
                                                  0, 1, m, x,
                                                  m, x, 0, n);
          
          auto bottom_right = Multiplier::multiply(input, H, // Okay since j can be considered constant in the maxing
                                                   0, 1, m, x,
                                                   m, x, n, x);
          
          // Skew the input to match precomputed table / LCS
          for (uint64_t i = 0; i < x; ++i)
            input(0, i) += i;
          auto top_left = Multiplier::multiply(input, H,
                                               0, 1, 0, m,
                                               0, m, 0, n);
          auto top_right = Multiplier::multiply(input, H,
                                                0, 1, 0, m,
                                                0, m, n, x);
          
          assert(top_left.getColumns() == n && top_right.getColumns() == x - n);
          assert(bottom_left.getColumns() == n && bottom_right.getColumns() == x - n);
          
          for (int64_t j = 0; j < n; ++j) {
            const int64_t top = top_left(0, j) - m;
            const int64_t bottom = bottom_left(0, j);
            
            O.push_back(max(top, bottom));
          }
          for (int64_t j = n; j < x; ++j) {
            const int64_t top = top_right(0, j - n) + n - j - m;
            const int64_t bottom = bottom_right(0, j - n) + n - j;
            
            O.push_back(max(top, bottom));
          }
        }
        
        return O;
      }
      
      static string name() {
        return "SimpleLCSDIST";
      }
      
    private:
      typedef MaxMultiply::Simple<Matrix<int64_t>> Multiplier;
      
      int64_t rows_;
      int64_t cols_;
      Matrix<int64_t> matrix_;
    };
    
    class PermutationDISTTable {
    public:
      PermutationDISTTable(int64_t rows, int64_t cols, vector<int64_t> row_specification)
        : rows_(rows), cols_(cols),
          row_map_(vector<int64_t>(row_specification.size(), 0)),
          column_map_(vector<int64_t>(row_specification.size(), 0))
      {
        for (int64_t i = 0; i < row_specification.size(); ++i) {
          assert(row_specification[i] < row_specification.size() && row_specification[i] >= 0);
          row_map_[i] = row_specification[i];
          column_map_[row_specification[i]] = i;
        }
      }
      
      int64_t rows() const { return rows_; }
      int64_t cols() const { return cols_; }
      
      void setRows(int64_t rows) { rows_ = rows; }
      void setCols(int64_t cols) { cols_ = cols; }
      
      int64_t getCol(int64_t row) const {
        return row_map_[row];
      }
      
      int64_t getRow(int64_t col) const {
        return column_map_[col];
      }
      
      bool operator()(int64_t i, int64_t j) const {
        return getCol(i) == j;
      }
      
      Matrix<int64_t> unfoldH() const {
        auto inverse_transform = [](const Matrix<int64_t>& A, int64_t rows) {
          Matrix<int64_t> result(A.getRows(), A.getColumns(), 0);
          for (int64_t i = 0; i < A.getRows(); ++i) {
            for (int64_t j = 0; j < A.getColumns(); ++j) {
              result(i, j) = j - (i - rows) - A(i, j);
            }
          }
          return result;
        };
        
        Matrix<int64_t> sigmaed(size() + 1, size() + 1, 0);
        for (int64_t i_ = 0; i_ <= size(); ++i_) {
          for (int64_t j_ = 0; j_ <= size(); ++j_) {
            
            int64_t sum = 0;
            for (int64_t i = i_; i < size(); ++i) {
              sum += (row_map_[i] < j_);
            }
            
            sigmaed(i_, j_) = sum;
          }
        }
        
        return inverse_transform(sigmaed, rows());
      }
      
      Matrix<int64_t> unfoldSigma() const {
        Matrix<int64_t> sigmaed(size() + 1, size() + 1, 0);
        for (int64_t i_ = 0; i_ <= size(); ++i_) {
          for (int64_t j_ = 0; j_ <= size(); ++j_) {
            
            int64_t sum = 0;
            for (int64_t i = i_; i < size(); ++i) {
              sum += (row_map_[i] < j_);
            }
            
            sigmaed(i_, j_) = sum;
          }
        }
        return sigmaed;
      }
      
      Matrix<int64_t> unfold() const {
        Matrix<int64_t> unfolded(size(), size(), 0);
        for (int64_t i = 0; i < size(); ++i)
          unfolded(i, row_map_[i]) = 1;
        return unfolded;
      }
      
      int64_t size() const {
        return row_map_.size();
      }
      
      unique_ptr<PermutationDISTTable> swapStrings() const {
        assert(rows() + cols() == size());
        
        vector<int64_t> transformed_rowspec(size(), 0);
        for (int64_t row = 0; row < size(); ++row) {
          assert(rows() + cols() - row - 1 >= 0 && rows() + cols() - row - 1 < size());
          transformed_rowspec[rows() + cols() - row - 1] = rows() + cols() - getCol(row) - 1;
        }
        return unique_ptr<PermutationDISTTable>(new PermutationDISTTable(cols(), rows(), transformed_rowspec));
      }
      
      vector<int64_t> apply(const vector<int64_t>& I) const {
        assert(this->size() + 1 == I.size());
        
        vector<int64_t> O; O.reserve(I.size());
        
        const int64_t x = I.size();
        const int64_t m = rows_;
        const int64_t n = cols_;
        
        {
          vector<int64_t> input = I;
          for (int64_t i = 0; i < m; ++i)
            input[i] += i - m;
          
          auto multiplied = maxmultiply(input, this);
          
          for (int64_t j = 0; j < x; ++j) {
            if (j >= n) {
              O.push_back(multiplied[j].value + n - j);
            } else {
              O.push_back(multiplied[j].value);
            }
          }
        }
        
        return O;
      }
      
      bool operator==(const PermutationDISTTable& other) const {
        return row_map_ == other.row_map_ && rows_ == other.rows_ && cols_ == other.cols_;
      }
      
      bool operator!=(const PermutationDISTTable& other) const {
        return !((*this) == other);
      }
      
      static unique_ptr<PermutationDISTTable> BaseCase(bool is_equal) {
        if (is_equal)
          return unique_ptr<PermutationDISTTable>(new PermutationDISTTable(1, 1, { 0, 1 }));
        else
          return unique_ptr<PermutationDISTTable>(new PermutationDISTTable(1, 1, { 1, 0 }));
      }
      
      static unique_ptr<PermutationDISTTable> concat(const PermutationDISTTable& A, const PermutationDISTTable& B) {
        vector<int64_t> concatted_row_map;
        concatted_row_map.reserve(A.size() + B.size());
        for (int64_t i = 0; i < A.size(); ++i)
          concatted_row_map.push_back(A.row_map_[i]);
        for (int64_t i = 0; i < B.size(); ++i)
          concatted_row_map.push_back(B.row_map_[i] + A.size());
        
        return unique_ptr<PermutationDISTTable>(new PermutationDISTTable(A.rows() + B.rows(), A.cols() + B.cols(), concatted_row_map));
      }
      
      static string name() {
        return "PermutationDISTTable";
      }
      
      struct MaxInfo {
        MaxInfo(int64_t value, int64_t pos) : value(value), pos(pos) { }
        int64_t value, pos;
        
        bool operator==(const MaxInfo& other) const {
          return value == other.value && pos == other.pos;
        }
      };
      
      static vector<MaxInfo> maxmultiply_slow(vector<int64_t> v, const PermutationDISTTable* dist) {
        assert(v.size() == dist->size() + 1);
        const int64_t x = v.size();
        
        vector<MaxInfo> result; result.reserve(x);
        
        const auto& H = dist->unfoldH();
        for (int64_t j = 0; j < x; ++j) {
          MaxInfo maximum(numeric_limits<int64_t>::min(), -1);
          for (int64_t i = 0; i < x; ++i) {
            const int64_t candidate = v[i] + H(i, j);
            if (maximum.value < candidate) {
              maximum.value = candidate;
              maximum.pos = i;
            }
          }
          
          result.push_back(maximum);
        }
        
        return result;
      }
      
      static vector<MaxInfo> maxmultiply(vector<int64_t> v, const PermutationDISTTable* dist) {
        assert(v.size() == dist->size() + 1);
        const int64_t x = v.size();
        
        vector<MaxInfo> result; result.reserve(x);
        {
          list<MaxInfo> candidates;
          vector<list<MaxInfo>::iterator> iterators(x, candidates.end());
          Utils::IntervalUnionFind interval(x);
          
          candidates.push_back(MaxInfo(v[x - 1] - (x - 1), x - 1));
          iterators[x - 1] = --(candidates.rbegin().base());
          for (int64_t i = x - 2; i >= 0; --i) {
            const int64_t value = v[i] - i;
            if (value > candidates.back().value) {
              candidates.back().value -= value;
              candidates.push_back(MaxInfo(value, i));
              iterators[i] = --(candidates.rbegin().base());
            } else {
              interval.Union(i, i + 1);
            }
          }
          
          auto locate = [&candidates] (int64_t k) {
            for (auto iter = candidates.rbegin(); iter != candidates.rend(); ++iter) {
              if (iter->pos > k) {
                if (iter == candidates.rbegin())
                  return candidates.end();
                
                return iter.base();
              }
            }
            return candidates.end();
          };
          
          // Base case
          // Max is in end of candidate list
          MaxInfo answer(candidates.back().value + dist->rows(), candidates.back().pos);
          result.push_back(answer);
          
          for (int64_t j = 1; j < x; ++j) {
            const int64_t k = dist->getRow(j - 1);
            const auto pos = interval.Find(k);
            
            // Find the rightmost element with pos <= k
            list<MaxInfo>::iterator t;
            if (pos == k) {
              t = iterators[pos];
            } else {
              assert(pos > k);
              if (candidates.size() == 0) {
                t = candidates.end();
              } else if (candidates.begin()->pos <= k) {
                t = candidates.begin(); // Start of candidate list
              } else if (iterators[pos] == --candidates.end()) {
                t = candidates.end(); // We are at the beginning of candidate list
              } else {
                t = iterators[pos]; ++t;
                assert(t->pos <= k);
              }
            }
            assert(t == locate(k));
            
            if (t != candidates.end()) {
              candidates.back().value--;
              
              if (t != candidates.begin()) {
                list<MaxInfo>::iterator t_next = t; --t_next;
                
                if (t_next != candidates.end()) {
                  t_next->value++;
                  if (t_next->value == 0) {
                    // Union intervals
                    iterators[t->pos] = candidates.end();
                    interval.Union(t->pos, t->pos + 1);
                    
                    if (t == --candidates.end()) {
                      (++candidates.rbegin())->value += candidates.back().value;
                      candidates.pop_back();
                    } else {
                      t_next->value = t->value;
                      candidates.erase(t);
                    }
                  }
                }
              }
            }
            
            /*
             auto print_t = [&](int64_t j) {
             cout << "t(" << j << "):";
             for (int64_t i = 0; i < P_sigma.getRows(); ++i) {
             if (j == -1)
             cout << "\t" << v[i] - i;
             else
             cout << "\t" << v[i] - i - P_sigma(i, j);
             }
             cout << endl << endl;
             };
             
             print_t(j);
             cout << "Candidates: " << endl;
             for (auto iter = candidates.rbegin(); iter != candidates.rend(); ++iter) {
             cout << "Value: " << iter->value << ", Pos: " << iter->pos << endl;
             }
             */
            
            // Max is in end of candidate list
            MaxInfo answer(candidates.back().value + j + dist->rows(), candidates.back().pos);
            result.push_back(answer);
          }
        }
        
#ifndef NDEBUG
        {
          const auto& H = dist->unfoldH();
          for (int64_t j = 0; j < x; ++j) {
            MaxInfo maximum(numeric_limits<int64_t>::min(), -1);
            for (int64_t i = 0; i < x; ++i) {
              const int64_t candidate = v[i] + H(i, j);
              if (maximum.value < candidate) {
                maximum.value = candidate;
                maximum.pos = i;
              }
            }
            
            assert(result[j].value == maximum.value);
            assert(v[result[j].pos] + H(result[j].pos, j) == maximum.value);
          }
        }
#endif
        
        return result;
      }
      
    private:
      int64_t rows_, cols_;
      vector<int64_t> row_map_, column_map_;
    };
    
    bool operator==(const PermutationDISTTable& permDIST, const SimpleLCSDISTTable& simpleDIST) {
      return permDIST.rows() == simpleDIST.rows() && permDIST.cols() == simpleDIST.cols() && permDIST.unfoldH() == simpleDIST.matrix();
    }
    
    bool operator!=(const PermutationDISTTable& permDIST, const SimpleLCSDISTTable& simpleDIST) {
      return !(permDIST == simpleDIST);
    }
    
    template <class DISTTable>
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
            Matrix<int64_t> matrix(m + n + 1, m + n + 1, numeric_limits<int64_t>::max());
            
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
                
                ::Simple::EditDistance calculator(a->associatedString.substr(a_start, a_stop - a_start),
                                                  b->associatedString.substr(b_start, b_stop - b_start));
                // cout << "Res: " << calculator.edit_distance() << endl;
                matrix(in, out) = calculator.calculate();
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
      vector<int64_t> apply(Production* A, Production* B, const vector<int64_t>& I) {
        assert(A->DISTTableIndex != -1 && B->DISTTableIndex != -1);
        assert(A->associatedString.size() + B->associatedString.size() + 1 == I.size());
        
        /*
         cout << "A: " << A->associatedString << endl << "B: " << B->associatedString << endl;
         cout << "Associated DIST:" << endl << repo_[A->DISTTableIndex][B->DISTTableIndex] << endl << endl;
         */
        
        return (*this)(A, B).apply(I);
      }
      
      const DISTTable& operator()(Production* a, Production* b) const {
        const int64_t i = a->DISTTableIndex;
        const int64_t j = b->DISTTableIndex;
        
        assert(i != -1 && j != -1 && i < repo_.size() && j < repo_[i].size());
        return repo_[i][j];
      }
      
      static string name() {
        return "SimpleDIST";
      }
      
    private:
      vector<vector<DISTTable>> repo_;
    };
    
    template <class DISTTable>
    class LCSDISTRepository {
    public:
      LCSDISTRepository(const StraightLineProgram& slpA,
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
            Matrix<int64_t> matrix(m + n + 1, m + n + 1, numeric_limits<int64_t>::max());
            
            string b_padded = string(m, '*') + b->associatedString + string(m, '*');
            // Fill up the DIST matrix
            for (int64_t in = 0; in <= m + n; in++) {
              for (int64_t out = 0; out <= m + n; out++) {
                // if (out < in - (m + n) / 2) {
                if (out < in - m) {
                  // matrix(in, out) = out - in + (m + n) / 2;
                  matrix(in ,out) = out - (in - m);
                } else {
                  const uint64_t iStart = in;
                  const uint64_t iStop = out + m;
                  ::Simple::LCS lcs_calculator(a->associatedString, b_padded.substr(iStart, iStop - iStart));
                  matrix(in, out) = lcs_calculator.calculate();
                }
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
      vector<int64_t> apply(Production* A, Production* B, const vector<int64_t>& I) {
        assert(A->DISTTableIndex != -1 && B->DISTTableIndex != -1);
        assert(A->associatedString.size() + B->associatedString.size() + 1 == I.size());
        
        return (*this)(A, B).apply(I);
      }
      
      const DISTTable& operator()(Production* a, Production* b) const {
        const int64_t i = a->DISTTableIndex;
        const int64_t j = b->DISTTableIndex;
        
        assert(i != -1 && j != -1 && i < repo_.size() && j < repo_[i].size());
        return repo_[i][j];
      }
      
      static string name() {
        return "PermutationDIST";
      }
      
    private:
      vector<vector<DISTTable>> repo_;
    };
    
    template <class DISTTable, class Merger>
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
            assert(!a->associatedString.empty());
            assert(!b->associatedString.empty());
            build(a, b);
            
            assert((*this)(a, b).size() == a->associatedString.size() + b->associatedString.size());
          }
        }
      }
      
      /**
       * Given two productions A, B and input I, return the output of applying the DIST table identified to A, B to I.
       */
      vector<int64_t> apply(Production* A, Production* B, const vector<int64_t>& I) {
        assert(A->DISTTableIndex != -1 && B->DISTTableIndex != -1);
        assert(A->associatedString.size() + B->associatedString.size() + 1 == I.size());
        
        return (*this)(A, B).apply(I);
      }
      
      const DISTTable& operator()(Production* a, Production* b) const {
        const int64_t i = a->DISTTableIndex;
        const int64_t j = b->DISTTableIndex;
        
        assert(i != -1 && j != -1 && i < repo_.getRows() && j < repo_.getColumns() && repo_(i, j) != nullptr);
        return *repo_(i, j);
      }
      
      static string name() {
        return "MergingDISTRepository(" + DISTTable::name() + ", " + Merger::name() + ")";
      }
      
    protected:
      shared_ptr<DISTTable> build(Terminal* a, Terminal* b)
      {
        // cout << "Term - Term" << endl;
        const int64_t is_equal = a->symbol() == b->symbol();
        auto result = shared_ptr<DISTTable>(DISTTable::BaseCase(is_equal));
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
          assert(a->type != KEY);
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
      
      int64_t a_dist_index_count_, b_dist_index_count_;
      Matrix<shared_ptr<DISTTable>> repo_;
    };
    
    class SimpleMerger {
    public:
      static shared_ptr<EditDistanceDISTTable> horizontalMerge(const EditDistanceDISTTable* left,
                                                               const EditDistanceDISTTable* right)
      {
        assert(left->rows() == right->rows());
        const int64_t rows = left->rows();
        const int64_t cols = left->cols() + right->cols();
        const int64_t N = rows + cols + 1;
        
        auto result = shared_ptr<EditDistanceDISTTable>(new EditDistanceDISTTable({ rows, cols, Matrix<int64_t>(N, N, 9) }));
        for (int64_t in = 0; in < N; ++in) {
          for (int64_t out = 0; out < N; ++out) {
            // Only in left table
            if (out <= left->cols()) {
              if (in < left->matrix().getRows())
                result->matrix()(in, out) = left->matrix()(in, out);
              else
                result->matrix()(in, out) = numeric_limits<int64_t>::max();
            } else if (in < left->matrix().getRows()) {
              // Spanning both tables
              int64_t left_in = in;
              int64_t right_out = out - left->cols();
              
              int64_t minimum = numeric_limits<int64_t>::max();
              for (int64_t right_in = right->rows(); right_in >= max(right_out - right->cols(), (int64_t)0); --right_in) {
                int64_t left_out = right_in + left->cols();
                
                int64_t left_path = left->matrix()(left_in, left_out);
                int64_t right_path = right->matrix()(right_in, right_out);
                if (left_path == numeric_limits<int64_t>::max() ||
                    right_path == numeric_limits<int64_t>::max())
                {
                  continue;
                }
                
                minimum = min(minimum, left_path + right_path);
              }
              
              result->matrix()(in, out) = minimum;
            } else {
              // Only in right table
              const int64_t right_in = in - left->matrix().getRows() + right->rows() + 1;
              const int64_t right_out = out - left->cols();
              result->matrix()(in, out) = right->matrix()(right_in, right_out);
            }
          }
        }
        
        return result;
      }
      
      static shared_ptr<EditDistanceDISTTable> verticalMerge(const EditDistanceDISTTable* top,
                                                             const EditDistanceDISTTable* bottom)
      {
        assert(top->cols() == bottom->cols());
        const int64_t rows = bottom->rows() + top->rows();
        const int64_t cols = bottom->cols();
        const int64_t N = rows + cols + 1;
        
        auto result = shared_ptr<EditDistanceDISTTable>(new EditDistanceDISTTable(rows, cols, Matrix<int64_t>(N, N, 9)));
        for (int64_t in = 0; in < N; ++in) {
          for (int64_t out = 0; out < N; ++out) {
            if (out >= bottom->matrix().getRows()) {
              // Only in top table
              if (in > bottom->rows()) {
                const int64_t top_in = in - bottom->rows();
                const int64_t top_out = out - bottom->rows();
                result->matrix()(in, out) = top->matrix()(top_in, top_out);
              } else {
                result->matrix()(in, out) = numeric_limits<int64_t>::max();
              }
            } else if (in > bottom->rows()) {
              // Spanning both
              int64_t minimum = numeric_limits<int64_t>::max();
              for (int64_t k = 0; k <= min(out, bottom->cols()); ++k) {
                const int64_t top_path = top->matrix()(in - bottom->rows(), k);
                const int64_t bottom_path = bottom->matrix()(bottom->rows() + k, out);
                if (top_path == numeric_limits<int64_t>::max() ||
                    bottom_path == numeric_limits<int64_t>::max())
                {
                  continue;
                }
                
                minimum = min(minimum, top_path + bottom_path);
              }
              result->matrix()(in, out) = minimum;
            } else {
              // Only in bottom table
              result->matrix()(in, out) = bottom->matrix()(in, out);
            }
          }
        }
        
        return result;
      }
      
      static string name() {
        return "SimpleMerger";
      }
    };
    
    class SimpleLCSDISTMerger {
    public:
      static shared_ptr<SimpleLCSDISTTable> horizontalMerge(const SimpleLCSDISTTable* left,
                                                            const SimpleLCSDISTTable* right)
      {
        auto transform = [](const SimpleLCSDISTTable* original) {
          Matrix<int64_t> transformed(original->matrix().getRows(), original->matrix().getRows(), 9);
          
          const int64_t m = original->rows();
          const int64_t n = original->cols();
          
          for (int64_t i = 0; i <= m; ++i) {
            for (int64_t i_ = i; i_ <= m; ++i_) {
              transformed(i + n, i_) = original->matrix()(m - i, m + n - i_) - m - i + i_;
            }
          }
          for (int64_t l = 1; l <= n; ++l) {
            for (int64_t i_ = 0; i_ <= m; ++i_) {
              assert(transformed(n - l, i_) == 9);
              transformed(n - l, i_) = original->matrix()(m + l, m + n - i_) - m + i_ + l;
            }
          }
          for (int64_t i = 0; i <= m; ++i) {
            for (int64_t l_ = 0; l_ < n; ++l_) {
              assert(transformed(i + n, m + n - l_) == 9);
              transformed(i + n, m + n - l_) = original->matrix()(m - i, l_) + n - l_ - i;
            }
          }
          for (int64_t l = 1; l <= n; ++l) {
            for (int64_t l_ = 0; l_ < n; ++l_) {
              assert(transformed(n - l, m + n - l_) == 9);
              transformed(n - l, m + n - l_) = original->matrix()(m + l, l_) + n + l - l_;
            }
          }
          for (int64_t i = n + 1; i <= m + n; ++i) {
            for (int64_t j = 0; j < i - n; ++j) {
              transformed(i, j) = j - i + n;
            }
          }
          
          return unique_ptr<SimpleLCSDISTTable>(new SimpleLCSDISTTable(n, m, move(transformed)));
        };
        
        auto merged = verticalMerge(transform(left).get(), transform(right).get());
        return shared_ptr<SimpleLCSDISTTable>(transform(merged.get()));
      }
      
      static shared_ptr<SimpleLCSDISTTable> verticalMerge(const SimpleLCSDISTTable* top,
                                                          const SimpleLCSDISTTable* bottom)
      {
        assert(top->cols() == bottom->cols());
        
        auto transform = [](const SimpleLCSDISTTable* table) {
          const int64_t x = table->matrix().getRows();
          Matrix<int64_t> result(x, x, 0);
          for (int64_t i = 0; i < x; ++i) {
            for (int64_t j = 0; j < x; ++j) {
              result(i, j) = j - (i - table->rows()) - table->matrix()(i, j);
            }
          }
          
          return result;
        };
        
        auto inverse_transform = [](const Matrix<int64_t>& A, int64_t rows) {
          Matrix<int64_t> result(A.getRows(), A.getColumns(), 0);
          for (int64_t i = 0; i < A.getRows(); ++i) {
            for (int64_t j = 0; j < A.getColumns(); ++j) {
              result(i, j) = j - (i - rows) - A(i, j);
            }
          }
          return result;
        };
        
        auto box = [](const Matrix<int64_t>& A) {
          Matrix<int64_t> result(A.getRows() - 1, A.getColumns() - 1, 0);
          for (int64_t i = 0; i < A.getRows() - 1; ++i) {
            for (int64_t j = 0; j < A.getColumns() - 1; ++j) {
              result(i, j) = A(i + 1, j) - A(i, j) - A(i + 1, j + 1) + A(i, j + 1);
            }
          }
          
          return result;
        };
        
        auto sigma = [](const Matrix<int64_t>& A) {
          Matrix<int64_t> result(A.getRows() + 1, A.getColumns() + 1, 0);
          for (int64_t i_ = 0; i_ < result.getRows(); ++i_) {
            for (int64_t j_ = 0; j_ < result.getColumns(); ++j_) {
              
              int64_t sum = 0;
              for (int64_t i = i_; i < A.getRows(); ++i) {
                for (int64_t j = 0; j < j_; ++j) {
                  sum += A(i, j);
                }
              }
              
              result(i_, j_) = sum;
            }
          }
          return result;
        };
        
        auto top_perm = box(transform(top));
        auto bottom_perm = box(transform(bottom));
        
        assert(top->rows() + top->cols() + 1 == top->matrix().getRows());
        
        const int64_t padded_size = top->rows() + bottom->rows() + top->cols();
        Matrix<int64_t> padded_top_perm(padded_size, padded_size, 0),
                        padded_bottom_perm(padded_size, padded_size, 0);
        
        const int64_t h_ = bottom->rows();
        
        for (int64_t i = 0; i < padded_size; ++i) {
          for (int64_t j = 0; j < padded_size; ++j) {
            // Padded top
            if (i < h_ && j < h_) {
              // TODO: Consider if this is correct! This is simply the identity matrix right now!
              // const int64_t N = h_;
              // const int64_t h = bottom->rows() % N;
              // const int64_t h = top->rows() % N;
              // padded_top_perm(i, j) = (j - i == h || j - i == h - N);
              padded_top_perm(i, j) = (j == i);
            } else if (i >= h_ && j >= h_) {
              padded_top_perm(i, j) = top_perm(i - h_, j - h_);
            }
            
            // Padded bottom
            if (i < bottom_perm.getRows() && j < bottom_perm.getColumns()) {
              padded_bottom_perm(i, j) = bottom_perm(i, j);
            } else if (i >= bottom_perm.getRows() && j >= bottom_perm.getColumns()) {
              // TODO: Consider if this is correct! This is simply the identity matrix right now!
              // const int64_t N = padded_size - bottom_perm.getRows();
              // const int64_t h = bottom->rows() % N;
              // const int64_t h = top->rows() % N;
              // cout << "h: " << h << endl;
              // padded_bottom_perm(i, j) = (j - i == h || j - i == h - N);
              padded_bottom_perm(i, j) = (j == i);
            }
          }
        }
        
        // cout << "Top:" << endl << padded_top_perm << endl << "Bottom:" << endl << padded_bottom_perm << endl;
        
        auto multiplied = MaxMultiply::Simple<Matrix<int64_t>>::minmultiply(sigma(padded_top_perm), sigma(padded_bottom_perm),
                                                                            0, padded_size + 1, 0, padded_size + 1,
                                                                            0, padded_size + 1, 0, padded_size + 1);
        
        const int64_t rows = bottom->rows() + top->rows();
        const int64_t cols = bottom->cols();
        auto result = shared_ptr<SimpleLCSDISTTable>(new SimpleLCSDISTTable(rows, cols, inverse_transform(multiplied, rows)));
        
        // cout << "Multiplied: " << endl << box(transform(result.get())) << endl;
        
        return result;
      }
      
      static string name() {
        return "SimpleLCSDISTMerger";
      }
    };
    
    class PermutationLCSMerger {
    public:
      static shared_ptr<PermutationDISTTable> horizontalMerge(const PermutationDISTTable* left,
                                                              const PermutationDISTTable* right)
      {
        auto merged = verticalMerge(left->swapStrings().get(), right->swapStrings().get());
        return shared_ptr<PermutationDISTTable>(merged->swapStrings());
      }
      
      static shared_ptr<PermutationDISTTable> verticalMerge(const PermutationDISTTable* top,
                                                            const PermutationDISTTable* bottom)
      {
        auto identity = [](int64_t size) -> unique_ptr<PermutationDISTTable> {
          vector<int64_t> id; id.reserve(size);
          for (int64_t i = 0; i < size; ++i)
            id.push_back(i);
          
          return unique_ptr<PermutationDISTTable>(new PermutationDISTTable(0, 0, id));
        };
        
        auto padded_top = PermutationDISTTable::concat(*identity(bottom->rows()), *top);
        auto padded_bottom = PermutationDISTTable::concat(*bottom, *identity(top->rows()));
        
        return shared_ptr<PermutationDISTTable>(minmultiply(padded_top.get(), padded_bottom.get()));
      }
      
      static string name() {
        return "PermutationLCSMerger";
      }
      
      static unique_ptr<PermutationDISTTable> minmultiply(const PermutationDISTTable* A,
                                                          const PermutationDISTTable* B,
                                                          const int64_t base_case_switch_size = 20)
      {
        assert(A->size() == B->size());
        
        vector<int64_t> id;
        for (int64_t i = 0; i < A->size(); ++i)
          id.push_back(i);
        
        auto base_case = [A,B](vector<int64_t> A_row_indices, vector<int64_t> A_col_indices,
                               vector<int64_t> B_row_indices, vector<int64_t> B_col_indices)
                    -> unique_ptr<PermutationDISTTable>
        {
          // auto time1 = chrono::high_resolution_clock::now();
          
          assert(A_row_indices.size() == B_row_indices.size());
          const uint64_t n = A_row_indices.size();
          
          vector<int64_t> A_row_specification; A_row_specification.reserve(n);
          for (auto row : A_row_indices) {
            for (int64_t col = 0; col < n; ++col) {
              if (A_col_indices[col] == A->getCol(row)) {
                A_row_specification.push_back(col);
                break;
              }
            }
          }
          assert(A_row_specification.size() == n);
          PermutationDISTTable A_unified(0, 0, A_row_specification);
          
          vector<int64_t> B_row_specification; B_row_specification.reserve(n);
          for (auto row : B_row_indices) {
            for (int64_t col = 0; col < n; ++col) {
              if (B_col_indices[col] == B->getCol(row)) {
                B_row_specification.push_back(col);
                break;
              }
            }
          }
          assert(B_row_specification.size() == n);
          PermutationDISTTable B_unified(0, 0, B_row_specification);
          
          auto multiplied = ::MaxMultiply::Simple<Matrix<int64_t>>::minmultiply(A_unified.unfoldSigma(), B_unified.unfoldSigma(),
                                                                                0, n + 1, 0, n + 1, 0, n + 1, 0, n + 1);
          
          // Fold the multiplied result
          vector<int64_t> result_row_spec; result_row_spec.reserve(n);
          // Matrix<int64_t> result(multiplied.getRows() - 1, multiplied.getColumns() - 1, 0);
          for (int64_t i = 0; i < multiplied.getRows() - 1; ++i) {
            for (int64_t j = 0; j < multiplied.getColumns() - 1; ++j) {
              if (multiplied(i + 1, j) - multiplied(i, j) - multiplied(i + 1, j + 1) + multiplied(i, j + 1) == 1) {
                result_row_spec.push_back(j);
                break;
              }
            }
          }
          assert(result_row_spec.size() == n);
          
          // auto time2 = chrono::high_resolution_clock::now();
          
          // cout << n << "\t" << chrono::duration_cast<chrono::nanoseconds>(time2 - time1).count() << endl;
          
          return unique_ptr<PermutationDISTTable>(new PermutationDISTTable(0, 0, result_row_spec));
        };
        
        function<unique_ptr<PermutationDISTTable>(vector<int64_t>,vector<int64_t>,vector<int64_t>,vector<int64_t>)> compute;
        compute = [&compute,&base_case,A,B,base_case_switch_size]
          (vector<int64_t> A_row_indices, vector<int64_t> A_col_indices,
           vector<int64_t> B_row_indices, vector<int64_t> B_col_indices)
                    -> unique_ptr<PermutationDISTTable>
        {
          assert(A_row_indices.size() == B_row_indices.size() &&
                 B_row_indices.size() == A_col_indices.size() &&
                 A_col_indices.size() == B_col_indices.size());
          const int64_t n = A_row_indices.size();

          // TODO: Consider larger base case
          if (n <= base_case_switch_size) {
            return move(base_case(A_row_indices, A_col_indices, B_row_indices, B_col_indices));
          } else {
            // auto time1 = chrono::high_resolution_clock::now();
            
            // Generate subproblems and remap indices
            const int64_t lo_subproblem_size = ceil_div<int64_t>(n, 2);
            const int64_t hi_subproblem_size = n - lo_subproblem_size;
            
            vector<int64_t> A_col_lo_indices(A_col_indices.begin(), A_col_indices.begin() + lo_subproblem_size);
            vector<int64_t> A_col_hi_indices(A_col_indices.begin() + lo_subproblem_size, A_col_indices.end());
            assert(A_col_hi_indices.size() == hi_subproblem_size);
            
            vector<int64_t> B_row_lo_indices(B_row_indices.begin(), B_row_indices.begin() + lo_subproblem_size);
            vector<int64_t> B_row_hi_indices(B_row_indices.begin() + lo_subproblem_size, B_row_indices.end());
            assert(B_row_hi_indices.size() == hi_subproblem_size);
            
            vector<int64_t> A_row_lo_large_to_small(n, -1), A_row_hi_large_to_small(n, -1);
            vector<int64_t> A_row_lo_small_to_large, A_row_hi_small_to_large;
            vector<int64_t> A_row_lo_indices, A_row_hi_indices;
            A_row_lo_indices.reserve(lo_subproblem_size); A_row_hi_indices.reserve(hi_subproblem_size);
            A_row_lo_small_to_large.reserve(lo_subproblem_size); A_row_hi_small_to_large.reserve(hi_subproblem_size);
            for (int64_t i = 0; i < n; ++i) {
              if (A->getCol(A_row_indices[i]) < A_col_hi_indices[0]) {
                A_row_lo_large_to_small[i] = A_row_lo_indices.size();
                A_row_lo_small_to_large.push_back(i);
                A_row_lo_indices.push_back(A_row_indices[i]);
              } else {
                A_row_hi_large_to_small[i] = A_row_hi_indices.size();
                A_row_hi_small_to_large.push_back(i);
                A_row_hi_indices.push_back(A_row_indices[i]);
              }
            }
            assert(A_row_lo_indices.size() == lo_subproblem_size && A_row_hi_indices.size() == hi_subproblem_size);
            
            vector<int64_t> B_col_lo_large_to_small(n, -1), B_col_hi_large_to_small(n, -1);
            vector<int64_t> B_col_lo_small_to_large, B_col_hi_small_to_large;
            vector<int64_t> B_col_lo_indices, B_col_hi_indices;
            B_col_lo_indices.reserve(lo_subproblem_size); B_col_hi_indices.reserve(hi_subproblem_size);
            B_col_lo_small_to_large.reserve(lo_subproblem_size); B_col_hi_small_to_large.reserve(hi_subproblem_size);
            for (int64_t j = 0; j < n; ++j) {
              if (B->getRow(B_col_indices[j]) < B_row_hi_indices[0]) {
                B_col_lo_large_to_small[j] = B_col_lo_indices.size();
                B_col_lo_small_to_large.push_back(j);
                B_col_lo_indices.push_back(B_col_indices[j]);
              } else {
                B_col_hi_large_to_small[j] = B_col_hi_indices.size();
                B_col_hi_small_to_large.push_back(j);
                B_col_hi_indices.push_back(B_col_indices[j]);
              }
            }
            assert(B_col_lo_indices.size() == lo_subproblem_size && B_col_hi_indices.size() == hi_subproblem_size);
            
            // auto time2 = chrono::high_resolution_clock::now();
            
            // Recursively solve subproblems
            const auto P_C_lo = compute(A_row_lo_indices, A_col_lo_indices, B_row_lo_indices, B_col_lo_indices);
            const auto P_C_hi = compute(A_row_hi_indices, A_col_hi_indices, B_row_hi_indices, B_col_hi_indices);
            
            // auto time3 = chrono::high_resolution_clock::now();
            
            // Invert index remap
            auto query_row_P_C_hi = [&](int64_t i) -> int64_t {
              if (A_row_hi_large_to_small[i] == -1) return -1;
              return B_col_hi_small_to_large[P_C_hi->getCol(A_row_hi_large_to_small[i])];
            };
            
            auto query_col_P_C_hi = [&](int64_t k) -> int64_t {
              if (B_col_hi_large_to_small[k] == -1) return -1;
              return A_row_hi_small_to_large[P_C_hi->getRow(B_col_hi_large_to_small[k])];
            };
            
            auto query_row_P_C_lo = [&](int64_t i) -> int64_t {
              if (A_row_lo_large_to_small[i] == -1) return -1;
              return B_col_lo_small_to_large[P_C_lo->getCol(A_row_lo_large_to_small[i])];
            };
            
            auto query_col_P_C_lo = [&](int64_t k) -> int64_t {
              if (B_col_lo_large_to_small[k] == -1) return -1;
              return A_row_lo_small_to_large[P_C_lo->getRow(B_col_lo_large_to_small[k])];
            };

#ifndef NDEBUG
            auto query_P_C_lo = [&](int64_t i, int64_t j) -> bool {
              if (A_row_lo_large_to_small[i] == -1 || B_col_lo_large_to_small[j] == -1) return false;
              return (*P_C_lo)(A_row_lo_large_to_small[i], B_col_lo_large_to_small[j]);
            };
            
            auto query_P_C_hi = [&](int64_t i, int64_t j) -> bool {
              if (A_row_hi_large_to_small[i] == -1 || B_col_hi_large_to_small[j] == -1) return false;
              return (*P_C_hi)(A_row_hi_large_to_small[i], B_col_hi_large_to_small[j]);
            };
            
            // Delta function (slow version)
            auto slow_delta = [&](int64_t i, int64_t k) -> int64_t {
              assert(i >= 0 && i <= n && k >= 0 && k <= n);
              int64_t result = 0;
              for (int64_t i_ = 0; i_ < i; ++i_) {
                for (int64_t k_ = 0; k_ < k; ++k_) {
                  result += query_P_C_hi(i_, k_);
                }
              }
              for (int64_t i_ = i; i_ < n; ++i_) {
                for (int64_t k_ = k; k_ < n; ++k_) {
                  result -= query_P_C_lo(i_, k_);
                }
              }
              
              return result;
            };
#endif
      
            /*
            cout << endl << endl << "Lo:" << endl;
            for (int64_t i = 0; i < n; ++i) {
              for (int64_t k = 0; k < n; ++k)
                cout << query_P_C_lo(i, k) << " ";
              cout << endl;
            }
            
            cout << "Hi:" << endl;
            for (int64_t i = 0; i < n; ++i) {
              for (int64_t k = 0; k < n; ++k)
                cout << query_P_C_hi(i, k) << " ";
              cout << endl;
            }
            
            cout << "Delta:" << endl;
            for (int64_t i = 0; i <= n; ++i) {
              for (int64_t k = 0; k <= n; ++k) {
                cout << slow_delta(i, k) << " ";
              }
              cout << endl;
            }
             */
            
            auto query_delta_up = [&](int64_t old_result, int64_t old_row, int64_t old_col) -> int64_t {
              int64_t delta = 0;
              if (old_row > 0) {
                const auto hi_row_query = query_row_P_C_hi(old_row - 1);
                if (hi_row_query != -1 && hi_row_query < old_col)
                  delta--;
                
                const auto lo_row_query = query_row_P_C_lo(old_row - 1);
                if (lo_row_query != -1 && lo_row_query >= old_col)
                  delta--;
              }
              return old_result + delta;
            };
            
            auto query_delta_right = [&](int64_t old_result, int64_t old_row, int64_t old_col) -> int64_t {
              int64_t delta = 0;
              if (old_col < n) {
                const auto hi_col_query = query_col_P_C_hi(old_col);
                if (hi_col_query != -1 && hi_col_query < old_row)
                  delta++;
                
                const auto lo_col_query = query_col_P_C_lo(old_col);
                if (lo_col_query != -1 && lo_col_query >= old_row)
                  delta++;
              }
              
              return old_result + delta;
            };
            
            // Compute seperating paths for negative and positive values in delta matrix
            const int64_t K = n + 1;
            struct R {
              int64_t r_high, r_low;
            };
            vector<R> r_(2 * K + 1, R());
            auto r = [&r_, K] (int64_t i) -> R& {
              assert(i + K >= 0 && i + K < 2 * K + 1);
              return r_[i + K];
            };
            
            // Base case
            r(-K).r_low = n;
            r(-K).r_high = n;
            
            // Inductive case
            int64_t high_witness_value = 0, low_witness_value = 0;
            for (int64_t d = -K; d < K; ++d) {
              /*
               * High path
               */
              // cout << "High pos: (" << double(r(d).r_high - d) / 2 << ", " << double(r(d).r_high + d) / 2 << ")" << endl;
              if (r(d).r_high - d - 1 < 0) { // We are at the top -> go right
                r(d + 1).r_high = r(d).r_high + 1;
              } else {
                // cout << "Witness pos: (" << double(r(d).r_high - d - 1) / 2 << ", " << double(r(d).r_high + d + 1) / 2 << ")" << endl;
                assert(high_witness_value == slow_delta((r(d).r_high - d - 1) / 2, (r(d).r_high + d + 1) / 2));
                
                if (high_witness_value < 0) { // Go right
                  r(d + 1).r_high = r(d).r_high + 1;
                  high_witness_value = query_delta_right(high_witness_value, (r(d).r_high - d - 1) / 2, (r(d).r_high + d + 1) / 2);
                } else { // Go up
                  r(d + 1).r_high = r(d).r_high - 1;
                  high_witness_value = query_delta_up(high_witness_value, (r(d).r_high - d - 1) / 2, (r(d).r_high + d + 1) / 2);
                }
              }
              
              /*
               * Low path
               */
              // cout << "Low pos: (" << double(r(d).r_low - d) / 2 << ", " << double(r(d).r_low + d) / 2 << ")" << endl;
              if (r(d).r_low + d + 1 >= 2 * K) { // We are at the right -> go up
                r(d + 1).r_low = r(d).r_low - 1;
              } else {
                // cout << "Witness pos: (" << double(r(d).r_low - d - 1) / 2 << ", " << double(r(d).r_low + d + 1) / 2 << ")" << endl;
                assert(low_witness_value == slow_delta((r(d).r_low - d - 1) / 2, (r(d).r_low + d + 1) / 2));
                
                if (low_witness_value > 0) { // Go up
                  r(d + 1).r_low = r(d).r_low - 1;
                  low_witness_value = query_delta_up(low_witness_value, (r(d).r_low - d - 1) / 2, (r(d).r_low + d + 1) / 2);
                } else { // Go right
                  r(d + 1).r_low = r(d).r_low + 1;
                  low_witness_value = query_delta_right(low_witness_value, (r(d).r_low - d - 1) / 2, (r(d).r_low + d + 1) / 2);
                }
              }
            }
            // cout << "High pos: (" << double(r(K).r_high - K) / 2 << ", " << double(r(K).r_high + K) / 2 << ")" << endl;
            assert(double(r(K).r_high - K) / 2 == -0.5 && double(r(K).r_high + K) / 2 == K - 0.5);
            // cout << "Low pos: (" << double(r(K).r_low - K) / 2 << ", " << double(r(K).r_low + K) / 2 << ")" << endl;
            assert(double(r(K).r_low - K) / 2 == -0.5 && double(r(K).r_low + K) / 2 == K - 0.5);
            
#ifndef NDEBUG
            for (int64_t i = 0; i < K; ++i) {
              for (int64_t k = 0; k < K; ++k) {
                {
                  if (slow_delta(i, k) >= 0)
                    assert(i + k >= r(k - i).r_high);
                  
                  if (i + k >= r(k - i).r_high)
                    assert(slow_delta(i, k) >= 0);
                }
                
                {
                  if (slow_delta(i, k) <= 0)
                    assert(i + k <= r(k - i).r_low);
                  
                  if (i + k <= r(k - i).r_low)
                    assert(slow_delta(i, k) <= 0);
                }
                
                {
                  if (slow_delta(i, k) < 0 && slow_delta(i + 1, k + 1) > 0) {
                    assert(r(k - i).r_low == r(k - i).r_high && r(k - i).r_high == i + k + 1);
                  }
                  
                  if (r(k - i).r_low == r(k - i).r_high && r(k - i).r_high == i + k + 1) {
                    assert(slow_delta(i, k) < 0 && slow_delta(i + 1, k + 1) > 0);
                  }
                }
              }
            }
#endif
            
            // auto time4 = chrono::high_resolution_clock::now();
            
            // Combine solution
            vector<int64_t> combined_row_description(n, -1);
            
            for (int64_t row = 0; row < lo_subproblem_size; ++row) {
              const int64_t i = A_row_lo_small_to_large[row];
              const int64_t k = B_col_lo_small_to_large[P_C_lo->getCol(row)];
              
              if (i + k + 2 <= r(k - i).r_low)
                combined_row_description[i] = k;
            }
            
            for (int64_t row = 0; row < hi_subproblem_size; ++row) {
              const int64_t i = A_row_hi_small_to_large[row];
              const int64_t k = B_col_hi_small_to_large[P_C_hi->getCol(row)];
              
              if (i + k >= r(k - i).r_high)
                combined_row_description[i] = k;
            }
            
            for (int64_t d = -n; d < n; ++d) {
              if (r(d).r_high == r(d).r_low) {
                const int64_t i = (r(d).r_high - d - 1) / 2;
                const int64_t k = (r(d).r_high + d - 1) / 2;
                assert(r(d).r_high == i + k + 1);
                
                combined_row_description[i] = k;
              }
            }
            
#ifndef NDEBUG
            for (int64_t i = 0; i < n; ++i) {
              for (int64_t k = 0; k < n; ++k) {
                if (slow_delta(i, k) < 0 && slow_delta(i + 1, k + 1) > 0) {
                  assert(combined_row_description[i] == k);
                }
              }
            }
#endif
            
            // auto time5 = chrono::high_resolution_clock::now();
            
            /*
            cout << "input\t" << n << "\t" << chrono::duration_cast<chrono::nanoseconds>(time2 - time1).count() << endl;
            cout << "traverse\t" << n << "\t" << chrono::duration_cast<chrono::nanoseconds>(time4 - time3).count() << endl;
            cout << "combine\t" << n << "\t" << chrono::duration_cast<chrono::nanoseconds>(time5 - time4).count() << endl;
            cout << "total\t" << n << "\t" << chrono::duration_cast<chrono::nanoseconds>(time5 - time1).count() << endl;
             */
            
            return unique_ptr<PermutationDISTTable>(new PermutationDISTTable(0, 0, combined_row_description));
          }
        };
        
        vector<int64_t> all_indices; all_indices.reserve(A->size());
        for (int64_t i = 0; i < A->size(); ++i)
          all_indices.push_back(i);
        auto result = compute(all_indices, all_indices, all_indices, all_indices);
        result->setRows(A->rows() + B->rows());
        result->setCols(A->cols());
        return result;
      }
    };
  }
}
