#pragma once

#include <vector>

#include "../utils/matrix.hpp"
#include "SLP.hpp"

namespace Compression {
  namespace DIST {
    using namespace std;
    using namespace SLP;
    
    class SimpleDISTRepository {
    public:
      SimpleDISTRepository() { }
      
      /**
       * Constructs DIST repository given blocks in both strings.
       */
      void build(const vector<Production*>& A, const vector<Production*>& B) {
        repo_.reserve(A.size());
        
        for (auto a : A) {
          vector<Matrix<uint64_t>> row;
          row.reserve(B.size());
          assert(a->DISTTableIndex == -1);
          a->DISTTableIndex = repo_.size();
          
          for (auto b : B) {
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
            
            row.push_back(move(matrix));
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
            
            if (repo_[A->DISTTableIndex][B->DISTTableIndex](i, j) == numeric_limits<uint64_t>::max()) continue;
            minimum = min(minimum, I[i] + repo_[A->DISTTableIndex][B->DISTTableIndex](i, j));
          }
          O.push_back(minimum);
        }
        
        return O;
      }
      
    private:
      vector<vector<Matrix<uint64_t>>> repo_;
    };
  }
}
