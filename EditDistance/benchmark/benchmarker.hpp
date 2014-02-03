#pragma once

#include <cstdint>
#include <string>
#include <functional>
#include <cmath>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <sstream>

using namespace std;
using namespace std::chrono;

namespace Benchmark {
  struct Measurement {
    Measurement() : time(0), l2_cache_hits(0), l2_cache_misses(0),
                    l3_cache_hits(0), l3_cache_misses(0), instructions(0)
    { }
    
    bool operator<(const Measurement& o) const {
      return time < o.time;
    }
    
    double time;
    uint64_t l2_cache_hits, l2_cache_misses,
             l3_cache_hits, l3_cache_misses,
             instructions;
  };
  
  struct Statistics {
    double RSD, mean;
  };
  
  template<class Return>
  Measurement time(function<Return()> fct) {
    Measurement measurement;
    
    auto start = high_resolution_clock::now();
    
    fct();
    
    auto duration = high_resolution_clock::now() - start;
    measurement.time = duration_cast<microseconds>(duration).count() / 1E6; // Time in seconds
    
    return measurement;
  }
  
  /**
   * NOTICE: Strings are NOT generated uniformly at random!!!
   */
  string generate_string(uint64_t length, vector<char> alphabet = { 'a', 'c', 'g', 't' })
  {
    stringstream ss;
    
    for (uint64_t i = 0; i < length; ++i)
      ss << alphabet[rand() % alphabet.size()]; // Not fair randomness, but its good enough for this testing!
    
    return ss.str();
  }
  
  template <class Implementation>
  void run_benchmark(uint16_t trials, const uint64_t x) {
    string a = generate_string(100);
    string b = generate_string(100);
    
    cout << "#Testing: " << Implementation::name() << endl;
    
    const auto& stages = Implementation::run({ "a", "b", 1 });
    for (auto& stage : stages) // Just run the implementation on some very short strings, since this is only to get names of stages!
      stage.second();
    
    vector<pair<vector<Measurement>, Statistics>> measurements(stages.size(), pair<vector<Measurement>, Statistics>());
    
    // Run the actual tests
    for (uint16_t trial = 0; trial < trials; ++trial) {
      auto run = Implementation::run({ a, b, x });
      for (uint16_t stage = 0; stage < run.size(); ++stage) {
        auto& fct = run[stage].second;
        
        Measurement measurement = time(fct);
        measurements[stage].first.push_back(measurement);
      }
    }
    
    // Compute interesting results from measurements
    for (uint16_t stage = 0; stage < stages.size(); ++stage) {
      sort(measurements[stage].first.begin(), measurements[stage].first.end());
      
      // Compute mean value
      double mean = 0;
      for (uint16_t trial = 0; trial < trials; trial++)
        mean += measurements[stage].first[trial].time;
      mean /= trials;
      
      // Compute %RSD
      double standard_deviation = 0;
      for (uint16_t trial = 0; trial < trials; trial++)
        standard_deviation += pow(measurements[stage].first[trial].time - mean, 2);
      standard_deviation = sqrt(standard_deviation / trials);
      
      measurements[stage].second.RSD = (standard_deviation / mean) * 100;
      measurements[stage].second.mean = mean;
    }
    
    const uint16_t iMin = 0;
    const uint16_t iMax = trials - 1;
    const uint16_t iLower = iMax / 4;
    const uint16_t iUpper = (3 * iMax) / 4;
    const uint16_t iMedian = trials / 2;
    
    // Print out the test results
    cout << left;
    cout << "#" << setw(14) << "Stage" << setw(10) << "Time [s]" << setw(60) << " " << endl; // << setw(15) << "Instructions" << endl;
    cout << "#" << setw(14) << " " << setw(10) << "min" << setw(10) << "lower" << setw(10) << "median" << setw(10) << "upper"
         << setw(10) << "max" << setw(10) << "mean" << setw(10) << "%RSD"
      // << setw(15) << "mean"
         << endl;
    for (uint16_t stage = 0; stage < stages.size(); ++stage) {
      cout << setw(15) << stages[stage].first
           << setw(10) << measurements[stage].first[iMin].time
           << setw(10) << measurements[stage].first[iLower].time
           << setw(10) << measurements[stage].first[iMedian].time
           << setw(10) << measurements[stage].first[iUpper].time
           << setw(10) << measurements[stage].first[iMax].time
           << setw(10) << measurements[stage].second.mean
           << setw(10) << measurements[stage].second.RSD
        // << setw(15) << measurements[stage].first[iMedian].instructions
           << endl;
    }
  }
};
