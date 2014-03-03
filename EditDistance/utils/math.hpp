#pragma once

namespace Utils {
  namespace Math {
    template <class T> T ceil_div(T a, T d) {
      return (a + (d - 1)) / d;
    }
    
    int msb(unsigned long long x) {
      if (x == 0) return -1;
      return sizeof(unsigned long long) * 8 - __builtin_clzll(x);
    }
  }
};
