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
    int64_t l2_cache_hits, l2_cache_misses,
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
  string generate_string(int64_t length, vector<char> alphabet = { 'a', 'c', 'g', 't' })
  {
    stringstream ss;
    
    for (int64_t i = 0; i < length; ++i)
      ss << alphabet[rand() % alphabet.size()]; // Not fair randomness, but its good enough for this testing!
    
    return ss.str();
  }
  
  string easy_compressible_string(int64_t length, vector<char> alphabet = { 'a', 'c', 'g', 't' }) {
    stringstream ss;
    for (int64_t i = 0; i < length; ++i)
      ss << alphabet[i % alphabet.size()];
    
    return ss.str();
  }
  
  string fib_string(uint64_t n) {
    vector<string> mem(n + 1, "");
    
    function<string(uint64_t)> fib;
    fib = [&fib,&mem](uint64_t n) -> string {
      if (!mem[n].empty())
        return mem[n];
      
      if (n == 0) return "a";
      else if (n == 1) return "ab";
      else {
        string res = fib(n - 1) + fib(n - 2);
        mem[n] = res;
        return res;
      }
    };
    return fib(n);
  };
  
  template <class Implementation>
  void run_benchmark(uint16_t trials, const double xfactor) {
    /*
    string a = generate_string(1000, { 'a' });
    string b = generate_string(1000, { 'a' });
     */
    
    string a = fib_string(23);
    string b = fib_string(23);
  
    cout << "#Testing: " << Implementation::name() << endl;
    
    const auto& stages = Implementation::run({ "a", "b", 1 });
    for (auto& stage : stages) // Just run the implementation on some very short strings, since this is only to get names of stages!
      stage.second();
    
    vector<pair<vector<Measurement>, Statistics>> measurements(stages.size(), pair<vector<Measurement>, Statistics>());
    
    // Run the actual tests
    for (uint16_t trial = 0; trial < trials; ++trial) {
      auto run = Implementation::run({ a, b, xfactor });
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
    double total_median_time = 0;
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
      
      total_median_time += measurements[stage].first[iMedian].time;
    }
    cout << setw(15) << "Total" << setw(20) << " " << setw(10) << total_median_time << endl;
    
    cout << endl;
  }
};
