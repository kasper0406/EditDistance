#pragma once

#include <string>
#include <utility>
#include <functional>
#include <vector>

using namespace std;

class Simple {
public:
  Simple(string a, string b) : a_(a), b_(b) { }
  
  static vector<pair<string, function<uint64_t()>>> run(tuple<string, string, uint64_t> input) {
    Simple* simple = new Simple();
    tie(simple->a_, simple->b_, ignore) = input;
    
    auto compute = [simple]() -> uint64_t {
      return simple->edit_distance();
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
  uint64_t edit_distance() {
    if (b_.size() < a_.size()) swap(a_, b_);
    if (a_.empty()) return b_.size();
    assert(a_.size() > 0 && b_.size() > 0);
    
    auto substitution = [](char a, char b) {
      return 1 - (a == b);
    };
    
    vector<uint64_t> prev(a_.size() + 1, 0), cur(a_.size() + 1, 0);
    
    // Base case
    for (uint64_t i = 0; i <= a_.size(); i++) {
      // cout << i << "\t";
      prev[i] = i;
    }
    // cout << endl;
    
    for (uint64_t j = 1; j <= b_.size(); j++) {
      for (uint64_t i = 0; i <= a_.size(); i++) {
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
    return "Simple";
  }
  
private:
  Simple() { }
  
  string a_, b_;
};
