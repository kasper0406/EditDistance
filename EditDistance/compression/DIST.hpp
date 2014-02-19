#pragma once

#include <vector>
#include <tuple>
#include <list>
#include <memory>

#include "../utils/matrix.hpp"
#include "SLP.hpp"
#include "maxmultiply.hpp"

namespace Compression {
  namespace DIST {
    typedef ::Compression::SLP::SLP StraightLineProgram;
    using namespace std;
    using namespace ::Compression::SLP;
    
    class EditDistanceDISTTable {
    public:
      EditDistanceDISTTable(int64_t rows, int64_t cols, Matrix<int64_t> matrix)
        : rows_(rows), cols_(cols), matrix_(move(matrix))
      { }
      
      int64_t rows() const { return rows_; }
      int64_t cols() const { return cols_; }
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
          row_map_[i] = row_specification[i];
          column_map_[row_specification[i]] = i;
        }
      }
      
      int64_t rows() const { return rows_; }
      int64_t cols() const { return cols_; }
      
      int64_t getCol(int64_t row) const {
        return row_map_[row];
      }
      
      int64_t getRow(int64_t col) const {
        return column_map_[col];
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
        vector<int64_t> transformed_rowspec(size(), 0);
        for (int64_t row = 0; row < size(); ++row)
          transformed_rowspec[rows() + cols() - row - 1] = rows() + cols() - getCol(row) - 1;
        return unique_ptr<PermutationDISTTable>(new PermutationDISTTable(cols(), rows(), transformed_rowspec));
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
            build(a, b);
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
              const int64_t N = h_;
              const int64_t h = bottom->rows() % N;
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
              const int64_t N = padded_size - bottom_perm.getRows();
              // const int64_t h = bottom->rows() % N;
              const int64_t h = top->rows() % N;
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
        auto identity = [](int64_t size) {
          vector<int64_t> id; id.reserve(size);
          for (int64_t i = 0; i < size; ++i)
            id.push_back(i);
          
          return PermutationDISTTable(0, 0, id);
        };
        
        auto padded_top = PermutationDISTTable::concat(identity(bottom->rows()), *top);
        auto padded_bottom = PermutationDISTTable::concat(*bottom, identity(top->rows()));
        
        return shared_ptr<PermutationDISTTable>(minmultiply(padded_top.get(), padded_bottom.get()));
      }
      
    private:
      static unique_ptr<PermutationDISTTable> minmultiply(const PermutationDISTTable* A,
                                                          const PermutationDISTTable* B)
      {
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
        
        auto box = [](const Matrix<int64_t>& A) {
          Matrix<int64_t> result(A.getRows() - 1, A.getColumns() - 1, 0);
          for (int64_t i = 0; i < A.getRows() - 1; ++i) {
            for (int64_t j = 0; j < A.getColumns() - 1; ++j) {
              result(i, j) = A(i + 1, j) - A(i, j) - A(i + 1, j + 1) + A(i, j + 1);
            }
          }
          
          return result;
        };
        
        auto result = box(MaxMultiply::Simple<Matrix<int64_t>>::minmultiply(sigma(A->unfold()), sigma(B->unfold()),
                                                                            0, A->size() + 1, 0, A->size() + 1,
                                                                            0, B->size() + 1, 0, B->size() + 1));
        
        vector<int64_t> compact_result(result.getRows(), 0);
        for (int64_t i = 0; i < result.getRows(); ++i) {
          for (int64_t j = 0; j < result.getColumns(); ++j) {
            if (result(i, j) == 1) {
              compact_result[i] = j;
              break;
            }
          }
        }
        
        return unique_ptr<PermutationDISTTable>(new PermutationDISTTable(A->rows() + B->rows(), A->cols(), compact_result));
      }
    };
  }
}
