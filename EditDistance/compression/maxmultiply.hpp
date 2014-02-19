#pragma once

#include <algorithm>
#include <string>
#include <memory>

namespace MaxMultiply {
  using namespace std;
  
  template <class Matrix>
  class Simple {
  public:
    /**
     * MaxMultiply submatrices of A and B specified by
     *  - rows: [rowStart, rowStop[
     *  - cols: [colStart, colStop[
     */
    static Matrix multiply(const Matrix& A, const Matrix& B,
                           int64_t rowStartA, int64_t rowStopA, int64_t colStartA, int64_t colStopA,
                           int64_t rowStartB, int64_t rowStopB, int64_t colStartB, int64_t colStopB) {
      const int64_t rowsA = rowStopA - rowStartA, colsA = colStopA - colStartA,
                    rowsB = rowStopB - rowStartB, colsB = colStopB - colStartB;
      
      assert(colsA == rowsB);
      Matrix result(rowsA, colsB, typename Matrix::Type(0));
      
      for (int64_t i = 0; i < rowsA; ++i) {
        for (int64_t j = 0; j < colsB; ++j) {
          typename Matrix::Type maximum = numeric_limits<typename Matrix::Type>::min();
          for (int64_t k = 0; k < colsA; ++k) {
            maximum = max(maximum, A(i + rowStartA, k + colStartA) + B(k + rowStartB, j + colStartB));
          }
          result(i, j) = maximum;
        }
      }
      
      return result;
    }
    
    /**
     * TODO: REFACTOR THIS!!!
     */
    static Matrix minmultiply(const Matrix& A, const Matrix& B,
                              int64_t rowStartA, int64_t rowStopA, int64_t colStartA, int64_t colStopA,
                              int64_t rowStartB, int64_t rowStopB, int64_t colStartB, int64_t colStopB) {
      const int64_t rowsA = rowStopA - rowStartA, colsA = colStopA - colStartA,
      rowsB = rowStopB - rowStartB, colsB = colStopB - colStartB;
      
      assert(colsA == rowsB);
      Matrix result(rowsA, colsB, typename Matrix::Type(0));
      
      for (int64_t i = 0; i < rowsA; ++i) {
        for (int64_t j = 0; j < colsB; ++j) {
          typename Matrix::Type minimum = numeric_limits<typename Matrix::Type>::max();
          for (int64_t k = 0; k < colsA; ++k) {
            minimum = min(minimum, A(i + rowStartA, k + colStartA) + B(k + rowStartB, j + colStartB));
          }
          result(i, j) = minimum;
        }
      }
      
      return result;
    }
    
    static string name() {
      return "SimpleMaxMultiplier";
    }
  };
};
