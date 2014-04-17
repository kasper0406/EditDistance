#pragma once

namespace Benchmark {
  struct Stats {
    Stats() : A_derivedLength(0), B_derivedLength(0), A_productions(0), B_productions(0), A_x(0), A_y(0), result(0) { }
    
    uint64_t A_derivedLength;
    uint64_t B_derivedLength;
    uint64_t A_productions;
    uint64_t B_productions;
    uint64_t A_x;
    uint64_t A_y;
    uint64_t result;
  };
};
