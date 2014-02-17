#pragma once

#include <string>
#include <utility>
#include <functional>
#include <vector>
#include <iostream>
#include <iomanip>

namespace Simple {
  using namespace std;
  
  class EditDistance {
  public:
    EditDistance(string a, string b) : a_(a), b_(b) { }
    
    static vector<pair<string, function<int64_t()>>> run(tuple<string, string, int64_t> input) {
      EditDistance* simple = new EditDistance();
      tie(simple->a_, simple->b_, ignore) = input;
      
      auto compute = [simple]() -> int64_t {
        return simple->calculate();
      };
      
      auto cleanup = [simple]() {
        delete simple;
        return 0;
      };
      
      return { { "computation", compute }, { "cleanup", cleanup } };
    }
    
    /**
     * TODO: Find a way of supplying a scoring function!
     */
    int64_t calculate() {
      if (b_.size() < a_.size()) swap(a_, b_);
      if (a_.empty()) return b_.size();
      assert(a_.size() > 0 && b_.size() > 0);
      
      auto substitution = [](char a, char b) {
        if (a == '*' || b == '*') return 0;
        
        return 1 - (a == b);
      };
      
      vector<int64_t> prev(a_.size() + 1, 0), cur(a_.size() + 1, 0);
      
      // Base case
      for (int64_t i = 0; i <= a_.size(); i++) {
        // cout << i << "\t";
        prev[i] = i;
      }
      // cout << endl;
      
      for (int64_t j = 1; j <= b_.size(); j++) {
        for (int64_t i = 0; i <= a_.size(); i++) {
          if (i == 0) {
            cur[i] = prev[i] + 1;
          } else {
            cur[i] = min(prev[i - 1] + substitution(a_[i - 1], b_[j - 1]),
                         min(prev[i] + 1, cur[i - 1] + 1));
          }
          
          // cout << cur[i] << "\t";
        }
        // cout << endl;
        
        swap(prev, cur);
      }
      
      return prev[a_.size()];
    }
    
    static string name() {
      return "Simple Edit-Distance";
    }
    
  private:
    EditDistance() { }
    
    string a_, b_;
  };
  
  class LCS {
  public:
    LCS(string a, string b) : a_(a), b_(b) { }
    
    int64_t calculate() {
      if (b_.size() < a_.size()) swap(a_, b_);
      if (a_.empty()) return 0;
      assert(a_.size() > 0 && b_.size() > 0);
      
      /*
      for (int64_t i = 0; i <= a_.size(); ++i)
        cout << left << setw(4) << 0;
      cout << endl;
       */
      
      vector<int64_t> prev(a_.size() + 1, 0), cur(a_.size() + 1, 0);
      for (int64_t j = 1; j <= b_.size(); j++) {
        for (int64_t i = 0; i <= a_.size(); i++) {
          if (i == 0) {
            cur[i] = prev[i];
          } else {
            cur[i] = max(prev[i - 1] + (a_[i - 1] == b_[j - 1] || a_[i - 1] == '*' || b_[j - 1] == '*'),
                         max(prev[i], cur[i - 1]));
          }
          
          // cout << left << setw(4) << cur[i];
        }
        // cout << endl;
        swap(prev, cur);
      }
      return prev[a_.size()];
    }
    
    static string name() {
      return "Simple LCS";
    }
    
  private:
    string a_, b_;
  };
  
  auto simple_outputting_lcs = [](string a_, string b_) -> int64_t {
    if (b_.size() < a_.size()) swap(a_, b_);
    if (a_.empty()) return 0;
    assert(a_.size() > 0 && b_.size() > 0);
    
    for (int64_t i = 0; i <= a_.size(); ++i)
      cout << left << setw(4) << 0;
    cout << endl;
    
    vector<int64_t> prev(a_.size() + 1, 0), cur(a_.size() + 1, 0);
    for (int64_t j = 1; j <= b_.size(); j++) {
      for (int64_t i = 0; i <= a_.size(); i++) {
        if (i == 0) {
          cur[i] = prev[i];
        } else {
          cur[i] = max(prev[i - 1] + (a_[i - 1] == b_[j - 1] || a_[i - 1] == '*' || b_[j - 1] == '*'),
                       max(prev[i], cur[i - 1]));
        }
        
        cout << left << setw(4) << cur[i];
      }
      cout << endl;
      swap(prev, cur);
    }
    return prev[a_.size()];
  };
}
