EditDistance
============

This repository provides the implementation for my Master's Thesis on accelerating the computation of Edit-Distance using compression.

The implemented algorithm is described in [CITE] and uses straight-line programs (SLPs) to encode the input strings.

The content of the repository include:
 - The Master's Thesis itself in the `report` folder.
 - The source code for the implemented algorithms in the `EditDistance` folder.
 - Raw measurements for all benchmarks discussed in the thesis. Located in the `Benchmarks` folder.
 - Plots and Python generators for all plots in the thesis. Also located in the `Benchmarks` folder.
 

### Running the program
All the implemented algorithms have been implemented in C++11 and has been compiled using the Clang 3.4 compiler. The implementation relies on the [SDSL lite](https://github.com/simongog/sdsl-lite) library and [Intel Performance Counter Monitor](https://software.intel.com/en-us/articles/intel-performance-counter-monitor-a-better-way-to-measure-cpu-utilization). The project uses `cmake` for building.

In order to build and run the implementation the following list of steps need to be completed:
 - Generate the make-file by running `cmake` in the `EditDistance` folder.
 - Build the implementation by running `make`.
 - Run the program by executing `./EditDistance`.

Depending on the computer configuration, it may be necessary to change the include and lib directories in the `CMakeLists.txt` file.

### Literature
TODO
