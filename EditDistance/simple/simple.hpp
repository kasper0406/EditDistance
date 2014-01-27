#pragma once

#include <string>
#include <utility>
#include <functional>
#include <vector>

using namespace std;

class Simple {
public:
  vector<pair<string, function<uint64_t()>>> run(pair<string, string> input) {
    a = input.first;
    b = input.second;
    
    auto compute = [&]() -> uint64_t {
      return edit_distance();
    };
    
    return { { "computation", compute } };
  }
  
  uint64_t edit_distance() {
    for (int i = 0; i < 100000000; ++i) continue;
    
    return 1;
  }
  
  static string name() {
    return "Simple";
  }
  
private:
  string a, b;
};
