#include <stdexcept>
#include <vector>

using namespace std;

template<class T>
class Matrix {
public:
  Matrix(size_t n, size_t m, const T& value)
    : n(n), m(m), elements(vector<T>(n * m, value))
  { }
  
  Matrix(size_t n, size_t m, vector<vector<T>> values)
    : Matrix(n, m, 0)
  {
    assert(values.size() == n);
    for (size_t i = 0; i < n; i++) {
      assert(values[i].size() == m);
      
      for (size_t j = 0; j < m; j++)
        (*this)(i, j) = values[i][j];
    }
  }
  
  Matrix(const Matrix<T>& other) {
    throw runtime_error("Do not copy a matrix!");
  }
  
  Matrix<T>& operator=(const Matrix<T>& other) {
    throw runtime_error("Do not copy a matrix!");
  }
  
  Matrix(const Matrix<T>&& other)
    : n(other.n), m(other.m), elements(move(other.elements))
  { }
  
  Matrix<T>& operator=(const Matrix<T>&& other)
  {
    n = other.n;
    m = other.m;
    elements = move(other.elements);
    return *this;
  }
  
  size_t getRows() const { return n; }
  size_t getColumns() const { return m; }
  
  inline typename vector<T>::const_reference operator()(size_t row, size_t column) const {
    assert(row < n && column < m);
    return elements[m * row + column];
  }
  
  inline typename vector<T>::reference operator()(size_t row, size_t column) {
    assert(row < n && column < m);
    return elements[m * row + column];
  }
  
private:
  size_t n, m;
  vector<T> elements;
};

template<class T>
ostream& operator<<(ostream& strm, const Matrix<T>& matrix)
{
  for (size_t i = 0; i < matrix.getRows(); i++) {
    for (size_t j = 0; j < matrix.getColumns(); j++) {
      if (matrix(i, j) == numeric_limits<T>::max())
        strm << "∞\t";
      else
        strm << matrix(i, j) << "\t";
    }
    strm << endl;
  }
  return strm;
}

template<class T>
bool operator==(const Matrix<T>& A, const Matrix<T>& B) {
  if (A.getRows() != B.getRows() || A.getColumns() != B.getColumns()) return false;
  
  bool all_equal = true;
  for (uint64_t row = 0; row < A.getRows(); ++row) {
    for (uint64_t col = 0; col < B.getColumns(); ++col)
      all_equal = all_equal && (A(row, col) == B(row, col));
  }
  
  return all_equal;
}

template<class T>
bool operator!=(const Matrix<T>& A, const Matrix<T>& B) {
  return !(A == B);
}
