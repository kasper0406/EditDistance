#pragma once

#include <cstdint>
#include <string>
#include <functional>
#include <cmath>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>

using namespace std;
using namespace std::chrono;

namespace Benchmark {
  struct Stats {
    Stats() : A_productions(0), B_productions(0), A_x(0), A_y(0) { }
    
    uint64_t A_productions;
    uint64_t B_productions;
    uint64_t A_x;
    uint64_t A_y;
  };
  
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
  }
  
  vector<pair<string,string>> read_fasta_from_stream(ifstream& stream) {
    vector<pair<string,string>> seqs;
    string seq, line, name;
    while (getline(stream, line)) {
      if (!line.empty() && line[0] == '>') {
        if (!seq.empty()) {
          seqs.push_back(make_pair(name, seq));
          seq = "";
        }
        name = line.substr(1);
      } else if (!line.empty() && line[0] == ';') {
        // Ignore the line...
      } else {
        for (unsigned int i = 0; i < line.length(); i++) {
          if (line[i] != ' ')
            seq += line[i];
        }
      }
    }
    seqs.push_back(make_pair(name, seq));
    return seqs;
  }
  
  vector<string> read_seqs_from_files(vector<string> files) {
    vector<string> seqs;
    for (string file : files) {
      ifstream stream(file);
      if (!stream.is_open())
        throw runtime_error("Could not find " + file);
      
      auto seqsInFile = read_fasta_from_stream(stream);
      for (auto seqdesc : seqsInFile)
        seqs.push_back(seqdesc.second);
      
      stream.close();
    }
    return seqs;
  }
  
  template <class Implementation>
  void run_benchmark(uint16_t trials, const double xfactor) {
    cout << "#Testing: " << Implementation::name() << endl;
    
    const auto& stages = Implementation::run({ "a", "b", 1 });
    for (auto& stage : stages) { // Just run the implementation on some very short strings, since this is only to get names of stages!
      stage.second();
    }
    
    // Output files
    vector<ofstream> outputs; outputs.reserve(stages.size() + 1);
    for (auto& stage : stages) {
      outputs.push_back(ofstream(Implementation::short_name() + "_" + stage.first + ".dat", ofstream::out));
    }
    outputs.push_back(ofstream(Implementation::short_name() + "_total.dat", ofstream::out));
    
    // Output headers
    for (auto& output : outputs) {
      output << left;
      output << "#" << setw(14) << "n" << setw(10) << "Time [s]" << setw(65) << " " << setw(20) << "productions" << setw(10) << "x" << endl; // << setw(15) << "Instructions" << endl;
      output << "#" << setw(14) << " " << setw(10) << "min" << setw(10) << "lower" << setw(10) << "median" << setw(10) << "upper"
             << setw(10) << "max" << setw(15) << "mean" << setw(10) << "%RSD"
             << setw(10) << "A" << setw(10) << "B" << setw(5) << "A" << setw(5) << "B"
             // << setw(15) << "mean"
             << endl;
    }
    
    const uint64_t maxN = 1000;
    for (uint64_t n = 10; n < maxN; n = 1.7 * n) {
      cout << "Testing n = " << n << endl;
      
      string a = read_seqs_from_files({ "genome1.fa" })[0].substr(0, n);
      string b = read_seqs_from_files({ "genome2.fa" })[0].substr(0, n);
      
      /*
       string a = generate_string(1000, { 'a' });
       string b = generate_string(1000, { 'a' });
       */
      
      /*
       string a = fib_string(23);
       string b = fib_string(23);
       */
      
      vector<pair<vector<Measurement>, Statistics>> measurements(stages.size(), pair<vector<Measurement>, Statistics>());
      
      Stats stats;
      // Run the actual tests
      for (uint16_t trial = 0; trial < trials; ++trial) {
        auto run = Implementation::run({ a, b, xfactor }, &stats);
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
      double total_median_time = 0;
      for (uint16_t stage = 0; stage < stages.size(); ++stage) {
        outputs[stage] << setw(15) << n
                       << setw(10) << measurements[stage].first[iMin].time
                       << setw(10) << measurements[stage].first[iLower].time
                       << setw(10) << measurements[stage].first[iMedian].time
                       << setw(10) << measurements[stage].first[iUpper].time
                       << setw(10) << measurements[stage].first[iMax].time
                       << setw(15) << measurements[stage].second.mean
                       << setw(10) << measurements[stage].second.RSD
                       << setw(10) << stats.A_productions
                       << setw(10) << stats.B_productions
                       << setw(5) << stats.A_x
                       << setw(5) << stats.A_y
                       // << setw(15) << measurements[stage].first[iMedian].instructions
                       << endl;
        
        total_median_time += measurements[stage].first[iMedian].time;
      }
      outputs[stages.size()] << setw(15) << n << setw(20) << " " << setw(10) << total_median_time
                             << setw(40) << " " << setw(10) << stats.A_productions << setw(10) << stats.B_productions
                             << setw(5) << stats.A_x << setw(5) << stats.A_y << endl;
    }
      
    for (auto& output : outputs)
      output.close();
  }
};
