#pragma once

#define IPCM

#include <cstdint>
#include <string>
#include <functional>
#include <cmath>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>

#include "utils.hpp"

#include "../compression/SLP.hpp"
#include "../compression/DIST.hpp"

#ifdef IPCM
#include "ipcm/cpucounters.h"
#endif

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
    
#ifdef IPCM
    SystemCounterState before_sstate = getSystemCounterState();
#endif
    
    auto start = high_resolution_clock::now();
    
    fct();
    
    auto duration = high_resolution_clock::now() - start;
    
#ifdef IPCM
    SystemCounterState after_sstate = getSystemCounterState();
    
    measurement.l2_cache_hits        = getL2CacheHits(before_sstate, after_sstate);
    measurement.l2_cache_misses      = getL2CacheMisses(before_sstate, after_sstate);
    measurement.l3_cache_hits        = getL3CacheHits(before_sstate, after_sstate);
    measurement.l3_cache_misses      = getL3CacheMisses(before_sstate, after_sstate);
    measurement.instructions         = getInstructionsRetired(before_sstate, after_sstate);
#endif
    
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
      
      if (n == 0) return "";
      if (n == 1) return "a";
      else if (n == 2) return "ab";
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
  
  class FibonacciInput {
  public:
    FibonacciInput(uint64_t to) : m_current(1), m_to(to) { }
    
    string name() const {
      return "fib";
    }
    
    bool hasNext() const {
      return m_current <= m_to;
    }
    
    string next() {
      return fib_string(++m_current);
    }
    
  private:
    uint64_t m_current, m_to;
  };
  
  class FastaInput {
  public:
    FastaInput(string file, uint64_t max_len = string::npos, double grow_factor = 1.7) : m_file(file), m_grow(grow_factor), m_current_len(1) {
      ifstream input(file);
      m_input = read_fasta_from_stream(input)[0].second.substr(0, max_len);
      input.close();
    }
    
    string name() const {
      return "fasta_" + m_file.substr(0, m_file.find("."));
    }
    
    bool hasNext() const {
      return m_current_len < m_input.size();
    }
    
    string next() {
      const string result = m_input.substr(0, m_current_len);
      m_current_len = max((uint64_t)ceil(m_grow * m_current_len), m_current_len + 1);
      return result;
    }
    
  private:
    string m_input;
    string m_file;
    double m_grow;
    uint64_t m_current_len;
  };
  
  class UniformRandomInput {
  public:
    UniformRandomInput(uint64_t max_len, double grow_factor = 1.7, vector<char> alphabet = { 'a', 'c', 'g', 't' }, uint64_t seed = 0xDEADBEEF)
      : m_max_len(max_len), m_current_len(1), m_grow_factor(grow_factor), m_alphabet(alphabet), m_generator(seed), m_distr(0, alphabet.size() - 1)
    { }
    
    string name() const {
      return "random";
    }
    
    bool hasNext() const {
      return m_current_len < m_max_len;
    }
    
    string next() {
      stringstream sstream;
      for (uint64_t i = 0; i < m_current_len; ++i)
        sstream << m_alphabet[m_distr(m_generator)];
      
      m_current_len = max((uint64_t)ceil(m_grow_factor * m_current_len), m_current_len + 1);
      return sstream.str();
    }
    
  private:
    uint64_t m_max_len, m_current_len;
    double m_grow_factor;
    vector<char> m_alphabet;
    
    mt19937 m_generator;
    uniform_int_distribution<int> m_distr;
  };
  
  template <class InputGenerator>
  void benchmark_compression(uint64_t trials, InputGenerator generator) {
    ofstream output("compression_" + generator.name() + ".dat", ofstream::out);
    
    output << left;
    output << setw(10) << "num" << setw(15) << "N" << setw(15) << "min" << setw(15) << "lower" << setw(15) << "median" << setw(15) << "upper"
                       << setw(15) << "max" << setw(15) << "compression";
    
#ifdef IPCM
    output << setw(15) << "L2_hits" << setw(15) << "L2_miss"
           << setw(15) << "L3_hits" << setw(15) << "L3_miss"
           << setw(15) << "instructions";
#endif
    output << endl;
    
    for (uint64_t iteration = 1; generator.hasNext(); ++iteration) {
      cout << "Running iteration: " << iteration << endl;
      
      auto input = generator.next();
      unique_ptr<Compression::SLP::SLP> slp;
      vector<Measurement> measurements;
      
      function<void()> compute = [&slp,input] () {
        slp = Compression::SLP::LZSLPBuilder::build(input);
      };
      for (uint64_t trial = 0; trial < trials; ++trial) {
        measurements.push_back(time(compute));
      }
      
      sort(measurements.begin(), measurements.end());
      const uint16_t iMin = 0;
      const uint16_t iMax = trials - 1;
      const uint16_t iLower = iMax / 4;
      const uint16_t iUpper = (3 * iMax) / 4;
      const uint16_t iMedian = trials / 2;
      
      output << setw(10) << iteration
             << setw(15) << input.size()
             << setw(15) << measurements[iMin].time
             << setw(15) << measurements[iLower].time
             << setw(15) << measurements[iMedian].time
             << setw(15) << measurements[iUpper].time
             << setw(15) << measurements[iMax].time
             << setw(15) << slp->compressionFactor();
      
#ifdef IPCM
      output << setw(15) << measurements[iMedian].l2_cache_hits
             << setw(15) << measurements[iMedian].l2_cache_misses
             << setw(15) << measurements[iMedian].l3_cache_hits
             << setw(15) << measurements[iMedian].l3_cache_misses
             << setw(15) << measurements[iMedian].instructions;
#endif
      
      output << endl;
    }
    
    output.close();
  }
  
  void benchmark_min_multiply(uint16_t trials) {
    using namespace Compression::DIST;
    
    ofstream output("min-multiply.dat", ofstream::out);
    
    output << left;
    output << setw(10) << "bcs" << setw(15) << "N" << setw(15) << "min" << setw(15) << "lower" << setw(15) << "median" << setw(15) << "upper"
           << setw(15) << "max";
    
#ifdef IPCM
    output << setw(15) << "L2_hits" << setw(15) << "L2_miss"
           << setw(15) << "L3_hits" << setw(15) << "L3_miss"
           << setw(15) << "instructions";
#endif
    output << endl;
    
    // vector<uint16_t> base_case_sizes = { 1, 10, 14, 17, 20, 24, 50 };
    vector<uint16_t> base_case_sizes = { 50 };
    // for (uint16_t bcs = 1; bcs < 50; ++bcs) {
    for (auto bcs : base_case_sizes) {
      cout << "Testing bcs = " << bcs << endl;
      
      for (uint64_t N = 1; N < 1000000; N = max(N + 1, (uint64_t)ceil(1.7 * N))) {
        cout << "Testing N = " << N << endl;
        
        vector<int64_t> A_rows(N, 0);
        vector<int64_t> B_rows(N, 0);
        
        for (int64_t i = 0; i < N; ++i) {
          A_rows[i] = i;
          B_rows[i] = i;
        }
        random_shuffle(A_rows.begin(), A_rows.end());
        random_shuffle(B_rows.begin(), B_rows.end());
        
        unique_ptr<PermutationDISTTable> A(new PermutationDISTTable(0, 0, A_rows));
        unique_ptr<PermutationDISTTable> B(new PermutationDISTTable(0, 0, B_rows));
        
        function<void()> compute = [&A,&B,bcs] () {
          auto result = PermutationLCSMerger::minmultiply(A.get(), B.get(), bcs);
        };
        vector<Measurement> measurements;
        for (uint64_t trial = 0; trial < trials; ++trial) {
          measurements.push_back(time(compute));
        }
        
        sort(measurements.begin(), measurements.end());
        const uint16_t iMin = 0;
        const uint16_t iMax = trials - 1;
        const uint16_t iLower = iMax / 4;
        const uint16_t iUpper = (3 * iMax) / 4;
        const uint16_t iMedian = trials / 2;
        
        output << setw(10) << bcs
               << setw(15) << N
               << setw(15) << measurements[iMin].time
               << setw(15) << measurements[iLower].time
               << setw(15) << measurements[iMedian].time
               << setw(15) << measurements[iUpper].time
               << setw(15) << measurements[iMax].time;
        
#ifdef IPCM
        output << setw(15) << measurements[iMedian].l2_cache_hits
               << setw(15) << measurements[iMedian].l2_cache_misses
               << setw(15) << measurements[iMedian].l3_cache_hits
               << setw(15) << measurements[iMedian].l3_cache_misses
               << setw(15) << measurements[iMedian].instructions;
#endif
        
        output << endl;
      }
    }
    
    output.close();
  }
  
  void benchmark_max_multiply(uint16_t trials) {
    using namespace Compression::DIST;
    
    ofstream output("max-multiply.dat", ofstream::out);
    
    output << left;
    output << setw(15) << "N" << setw(15) << "min" << setw(15) << "lower" << setw(15) << "median" << setw(15) << "upper"
           << setw(15) << "max";
    
#ifdef IPCM
    output << setw(15) << "L2_hits" << setw(15) << "L2_miss"
           << setw(15) << "L3_hits" << setw(15) << "L3_miss"
           << setw(15) << "instructions";
#endif
    output << endl;
    
    for (uint64_t N = 1; N < 1000000; N = max(N + 1, (uint64_t)ceil(1.7 * N))) {
      cout << "Testing N = " << N << endl;
      
      vector<int64_t> vec(N, 0);
      for (int64_t i = 0; i < N; ++i) vec[i] = i;
      random_shuffle(vec.begin(), vec.end());
      
      vector<int64_t> perm_rows(N, 0);
      for (int64_t i = 0; i < N; ++i) perm_rows[i] = i;
      random_shuffle(perm_rows.begin(), perm_rows.end());
      unique_ptr<PermutationDISTTable> perm(new PermutationDISTTable(0, 0, perm_rows));
      
      function<void()> compute = [&vec,&perm] () {
        auto result = PermutationDISTTable::maxmultiply(vec, perm.get());
      };
      vector<Measurement> measurements;
      for (uint64_t trial = 0; trial < trials; ++trial) {
        measurements.push_back(time(compute));
      }
      
      sort(measurements.begin(), measurements.end());
      const uint16_t iMin = 0;
      const uint16_t iMax = trials - 1;
      const uint16_t iLower = iMax / 4;
      const uint16_t iUpper = (3 * iMax) / 4;
      const uint16_t iMedian = trials / 2;
      
      output << setw(15) << N
             << setw(15) << measurements[iMin].time
             << setw(15) << measurements[iLower].time
             << setw(15) << measurements[iMedian].time
             << setw(15) << measurements[iUpper].time
             << setw(15) << measurements[iMax].time;
      
#ifdef IPCM
      output << setw(15) << measurements[iMedian].l2_cache_hits
             << setw(15) << measurements[iMedian].l2_cache_misses
             << setw(15) << measurements[iMedian].l3_cache_hits
             << setw(15) << measurements[iMedian].l3_cache_misses
             << setw(15) << measurements[iMedian].instructions;
#endif
      
      output << endl;
    }
    
    output.close();
  }
  
  void benchmark_slow_max_multiply(uint16_t trials) {
    using namespace Compression::DIST;
    
    ofstream output("slow-max-multiply.dat", ofstream::out);
    
    output << left;
    output << setw(15) << "N" << setw(15) << "min" << setw(15) << "lower" << setw(15) << "median" << setw(15) << "upper"
    << setw(15) << "max";
    
#ifdef IPCM
    output << setw(15) << "L2_hits" << setw(15) << "L2_miss"
    << setw(15) << "L3_hits" << setw(15) << "L3_miss"
    << setw(15) << "instructions";
#endif
    output << endl;
    
    for (uint64_t N = 1; N < 1000000; N = max(N + 1, (uint64_t)ceil(1.7 * N))) {
      cout << "Testing N = " << N << endl;
      
      vector<int64_t> vec(N, 0);
      for (int64_t i = 0; i < N; ++i) vec[i] = i;
      random_shuffle(vec.begin(), vec.end());
      
      vector<int64_t> perm_rows(N, 0);
      for (int64_t i = 0; i < N; ++i) perm_rows[i] = i;
      random_shuffle(perm_rows.begin(), perm_rows.end());
      unique_ptr<PermutationDISTTable> perm(new PermutationDISTTable(0, 0, perm_rows));
      
      function<void()> compute = [&vec,&perm] () {
        auto result = PermutationDISTTable::maxmultiply_slow(vec, perm.get());
      };
      vector<Measurement> measurements;
      for (uint64_t trial = 0; trial < trials; ++trial) {
        measurements.push_back(time(compute));
      }
      
      sort(measurements.begin(), measurements.end());
      const uint16_t iMin = 0;
      const uint16_t iMax = trials - 1;
      const uint16_t iLower = iMax / 4;
      const uint16_t iUpper = (3 * iMax) / 4;
      const uint16_t iMedian = trials / 2;
      
      output << setw(15) << N
      << setw(15) << measurements[iMin].time
      << setw(15) << measurements[iLower].time
      << setw(15) << measurements[iMedian].time
      << setw(15) << measurements[iUpper].time
      << setw(15) << measurements[iMax].time;
      
#ifdef IPCM
      output << setw(15) << measurements[iMedian].l2_cache_hits
      << setw(15) << measurements[iMedian].l2_cache_misses
      << setw(15) << measurements[iMedian].l3_cache_hits
      << setw(15) << measurements[iMedian].l3_cache_misses
      << setw(15) << measurements[iMedian].instructions;
#endif
      
      output << endl;
    }
    
    output.close();
  }
  
  template <class Implementation, class InputGenerator>
  void run_benchmark(uint16_t trials, const double xfactor, InputGenerator generator1, InputGenerator generator2) {
    cout << "#Testing: " << Implementation::name() << endl;
    
    const auto& stages = Implementation::run({ "a", "b", 1 });
    for (auto& stage : stages) { // Just run the implementation on some very short strings, since this is only to get names of stages!
      stage.second();
    }
    
    // Output files
    vector<ofstream> outputs; outputs.reserve(stages.size() + 1);
    for (auto& stage : stages) {
      outputs.push_back(ofstream(Implementation::short_name() + "_" + generator1.name() + "_" + generator2.name() + "_" + stage.first + ".dat", ofstream::out));
    }
    outputs.push_back(ofstream(Implementation::short_name() + "_" + generator1.name() + "_" + generator2.name() + "_total.dat", ofstream::out));
    
    // Output headers
    for (auto& output : outputs) {
      /*
      output << left;
      output << "#" << setw(14) << "n" << setw(10) << "Time [s]" << setw(65) << " " << setw(20) << "productions" << setw(10) << "x" << endl; // << setw(15) << "Instructions" << endl;
      output << "#" << setw(14) << " " << setw(10) << "min" << setw(10) << "lower" << setw(10) << "median" << setw(10) << "upper"
             << setw(10) << "max" << setw(15) << "mean" << setw(10) << "%RSD"
             << setw(10) << "A" << setw(10) << "B" << setw(5) << "A" << setw(5) << "B"
             // << setw(15) << "mean"
             << endl;
       */
      
      output << left;
      output << setw(15) << "iteration" << setw(15) << "A_len" << setw(15) << "B_len" << setw(15) << "min" << setw(15) << "lower" << setw(15) << "median" << setw(15) << "upper"
             << setw(15) << "max" << setw(15) << "mean" << setw(10) << "%RSD"
             << setw(10) << "A_prod" << setw(10) << "B_prod" << setw(10) << "A_x" << setw(10) << "B_x";
      
#ifdef IPCM
      output << setw(15) << "L2_hits" << setw(15) << "L2_miss"
             << setw(15) << "L3_hits" << setw(15) << "L3_miss"
             << setw(15) << "instructions";
#endif
      
      output << endl;
    }
    
    string a = "";
    string b = "";
    
    for (uint64_t iteration = 1; generator1.hasNext() || generator2.hasNext(); ++iteration) {
      cout << "Running iteration: " << iteration << endl;
      
      if (generator1.hasNext())
        a = generator1.next();
      if (generator2.hasNext())
        b = generator2.next();
      
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
      
      cout << "Result: " << stats.result << endl;
      
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
      uint64_t l2_misses = 0, l2_hits = 0, l3_misses = 0, l3_hits = 0, instructions = 0;
      for (uint16_t stage = 0; stage < stages.size(); ++stage) {
        outputs[stage] << setw(15) << iteration
                       << setw(15) << a.length()
                       << setw(15) << b.length()
                       << setw(15) << measurements[stage].first[iMin].time
                       << setw(15) << measurements[stage].first[iLower].time
                       << setw(15) << measurements[stage].first[iMedian].time
                       << setw(15) << measurements[stage].first[iUpper].time
                       << setw(15) << measurements[stage].first[iMax].time
                       << setw(15) << measurements[stage].second.mean
                       << setw(10) << measurements[stage].second.RSD
                       << setw(10) << stats.A_productions
                       << setw(10) << stats.B_productions
                       << setw(10) << stats.A_x
                       << setw(10) << stats.A_y;
                       // << setw(15) << measurements[stage].first[iMedian].instructions
        
#ifdef IPCM
        outputs[stage] << setw(15) << measurements[stage].first[iMedian].l2_cache_hits
                       << setw(15) << measurements[stage].first[iMedian].l2_cache_misses
                       << setw(15) << measurements[stage].first[iMedian].l3_cache_hits
                       << setw(15) << measurements[stage].first[iMedian].l3_cache_misses
                       << setw(15) << measurements[stage].first[iMedian].instructions;
#endif
        
        outputs[stage] << endl;
        
        total_median_time += measurements[stage].first[iMedian].time;
        l2_misses += measurements[stage].first[iMedian].l2_cache_misses;
        l2_hits += measurements[stage].first[iMedian].l2_cache_hits;
        l3_misses += measurements[stage].first[iMedian].l3_cache_misses;
        l3_hits += measurements[stage].first[iMedian].l3_cache_hits;
        instructions += measurements[stage].first[iMedian].instructions;
      }
      outputs[stages.size()] << setw(15) << iteration
                             << setw(15) << a.length()
                             << setw(15) << b.length()
                             << setw(15) << 0
                             << setw(15) << 0
                             << setw(15) << total_median_time
                             << setw(15) << 0
                             << setw(15) << 0
                             << setw(15) << 0
                             << setw(10) << 0
                             << setw(10) << stats.A_productions
                             << setw(10) << stats.B_productions
                             << setw(10) << stats.A_x
                             << setw(10) << stats.A_y;
                             // << setw(15) << measurements[stage].first[iMedian].instructions
      
#ifdef IPCM
      outputs[stages.size()] << setw(15) << l2_hits
                             << setw(15) << l2_misses
                             << setw(15) << l3_hits
                             << setw(15) << l3_misses
                             << setw(15) << instructions;
#endif
      
      outputs[stages.size()] << endl;
    }
    
    for (auto& output : outputs)
      output.close();
  }
};
