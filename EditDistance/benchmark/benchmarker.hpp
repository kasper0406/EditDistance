#pragma once

#include <cstdint>
#include <string>
#include <functional>

using namespace std;

namespace Benchmark {
  struct Measurement {
    Measurement() : time(0), l2_cache_hits(0), l2_cache_misses(0),
                    l3_cache_hits(0), l3_cache_misses(0), instructions(0),
                    weight(1)
    { }
    
    bool operator<(const Measurement& o) const {
      return time < o.time;
    }
    
    double time;
    uint64_t l2_cache_hits, l2_cache_misses,
             l3_cache_hits, l3_cache_misses,
             instructions;
    uint16_t weight;
  };
  
  struct Statistics {
    double RSD, mean;
  };
  
  template<class Return>
  Measurement time(function<Return()> fct) {
    Measurement Measurement;
    Measurement.time = 1;
    
    fct();
    
    return Measurement;
  }
  
  template <class Implementation>
  void run_benchmark(uint16_t trials) {
    string a = "aaaaaaaa";
    string b = "aaaaaaaa";
    
    cout << "Testing: " << Implementation::name() << endl;
    
    const auto& stages = Implementation().run({ a, b });
    vector<pair<vector<Measurement>, Statistics>> measurements(stages.size(), pair<vector<Measurement>, Statistics>());
    
    // Run the actual tests
    for (uint16_t trial = 0; trial < trials; ++trial) {
      auto run = Implementation().run({ a, b });
      for (uint16_t stage = 0; stage < stages.size(); ++stage) {
        auto& fct = run[stage].second;
        
        Measurement measurement = time(fct);
        measurements[stage].first.push_back(measurement);
      }
    }
    
    // Compute interesting results from measurements
    for (uint16_t stage = 0; stage < stages.size(); ++stage) {
      sort(measurements[stage].first.begin(), measurements[stage].first.end());
    }
    
    const uint16_t iMin = 0;
    const uint16_t iMax = trials - 1;
    const uint16_t iLower = iMax / 4;
    const uint16_t iUpper = (3 * iMax) / 4;
    const uint16_t iMedian = trials / 2;
    
    // Print out the test results
  }
};
