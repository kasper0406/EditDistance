#pragma once

#include "SLP.hpp"

namespace DIST {
  using namespace SLP;
  
  class SimpleDIST {
  public:
    SimpleDIST(const vector<pair<string, Production*>>& blocks) {
      // TODO: Do the actual computation
      
      uint64_t count = 0;
      for (auto& block : blocks) {
        get<1>(block)->DISTTableIndex = count++;
      }
    }
    
  private:
  };
}
