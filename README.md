EditDistance
============

This repository provides the implementation for my Master's Thesis on accelerating the computation of Edit-Distance using compression.

The implemented algorithm is described in [9] and uses straight-line programs (SLPs) to encode the input strings.

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

### Structure of code
The following is a list with a short description of the most interesting code files:

| File | Description |
| --- | --- |
| main.cpp | Driver code for running tests and computations. |
| test/tests.hpp | Unit-tests for the implemented algorithms. |
| utils/unionfind.hpp | Implementation of the interval union-find algorithm |
| benchmark/benchmarker.hpp | Code for running benchmarks and measuring the running time / CPU intrinsics |
| simple/simple.hpp | Implementation of the simple algorithm for computing Edit-Distance and LCS. |
| compression/SLP.hpp | Code for representing, building and doing operations on straight-line programs (SLPs). |
| compression/DIST.hpp | Code for building and using the DIST repository. |
| compression/slpalign.hpp | Code doing the actual Edit-Distance calculation by utilizing the DIST repository and SLPs |

### Literature
1. Nielsen, Dehouck, and Raza (2013) Burrows (1994) Chen and Chao (2010) Alexandre Tiskin (2013) Hermelin et al. (2010) Gawrychowski (2012) Bioinformatics (2013) Itai (2006) Espeholt (2013) Ohlebusch, Fischer, and Gog (2010) Rytter (2003) S. e. a. Gog (2014) Sakamoto (2005) Alexander Tiskin (2010) Gusfield (1997)

2. Bioinformatics, UCSC Genome. 2013. “GRCh38 Genome Reference Consortium Human Reference 38.” <http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/>. December. <http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/>.

3. Burrows, David J., Michael; Wheeler. 1994. “A Block Sorting Lossless Data Compression Algorithm.” *Digital Equipment Corporation, Technical Report 124*. <http://www.hpl.hp.com/techreports/Compaq-DEC/SRC-RR-124.pdf>.

4. Chen, Kuan-Yu, and Kun-Mao Chao. 2010. “A Fully Compressed Algorithm for Computing the Edit Distance of Run-Length Encoded Strings.” In *Proceedings of the 18th Annual European Conference on Algorithms: Part I*, 415–26. ESA’10. Berlin, Heidelberg: Springer-Verlag. <http://dl.acm.org/citation.cfm?id=1888935.1888984>.

5. Espeholt, Lasse. 2013. “Exploring the Practicality of the Four Russians Method for Sequence Alignment.” <https://docs.google.com/file/d/0B1ZO_s5g90zHVEtiNXpMV0ZUUEE>. <https://docs.google.com/file/d/0B1ZO_s5g90zHVEtiNXpMV0ZUUEE>.

6. Gawrychowski, Paweł. 2012. “Faster Algorithm for Computing the Edit Distance Between SLP-Compressed Strings.” In *Proceedings of the 19th International Conference on String Processing and Information Retrieval*, 229–36. SPIRE’12. Berlin, Heidelberg: Springer-Verlag. doi:[10.1007/978-3-642-34109-0\_24](http://dx.doi.org/10.1007/978-3-642-34109-0_24). <http://dx.doi.org/10.1007/978-3-642-34109-0_24>.

7. Gog, Simon et. al. 2014. “Succinct Data Structure Library (SDSL) 2.0.” <https://github.com/simongog/sdsl-lite>. March. <https://github.com/simongog/sdsl-lite>.

8. Gusfield, D. 1997. *Algorithms on Strings, Trees and Sequences: Computer Science and Computational Biology*. Cambridge University Press. <http://books.google.dk/books?id=Ofw5w1yuD8kC>.

9. Hermelin, Danny, Gad M. Landau, Shir Landau, and Oren Weimann. 2010. “Unified Compression-Based Acceleration of Edit-Distance Computation.” *CoRR* abs/1004.1194.

10. Itai, Alon. 2006. “Linear Time Restricted Union/Find.”

11. Nielsen, Kasper, Mathieu Dehouck, and Sarfraz Raza. 2013. “Advanced Algorithms Project 1 – Fibonacci Heaps.” <https://github.com/kasper0406/aa13/tree/master/Project1>. <https://github.com/kasper0406/aa13/tree/master/Project1>.

12. Ohlebusch, Enno, Johannes Fischer, and Simon Gog. 2010. “CST++.” In *String Processing and Information Retrieval - 17th International Symposium (SPIRE 2010)*, 322–33. <http://dx.doi.org/10.1007/978-3-642-16321-0_34>.

13. Rytter, Wojciech. 2003. “Application of Lempel–Ziv Factorization to the Approximation of Grammar-Based Compression.” *Theoretical Computer Science* 302 (1–3): 211–22. doi:[http://dx.doi.org/10.1016/S0304-3975(02)00777-6](http://dx.doi.org/http://dx.doi.org/10.1016/S0304-3975(02)00777-6). <http://www.sciencedirect.com/science/article/pii/S0304397502007776>.

14. Sakamoto, Hiroshi. 2005. “A Fully Linear-Time Approximation Algorithm for Grammar-Based Compression.” *Journal of Discrete Algorithms* 3 (2–4): 416–30. doi:[http://dx.doi.org/10.1016/j.jda.2004.08.016](http://dx.doi.org/http://dx.doi.org/10.1016/j.jda.2004.08.016). <http://www.sciencedirect.com/science/article/pii/S1570866704000632>.

15. Tiskin, Alexander. 2010. “Fast Distance Multiplication of Unit-Monge Matrices.” In *Proceedings of the Twenty-First Annual ACM-SIAM Symposium on Discrete Algorithms*, 1287–96. SODA ’10. Philadelphia, PA, USA: Society for Industrial; Applied Mathematics. <http://dl.acm.org/citation.cfm?id=1873601.1873704>.

16. Tiskin, Alexandre. 2013. “Semi-Local String Comparison: algorithmic Techniques and Applications.” *CoRR* abs/0707.3619v21. <http://arxiv.org/pdf/0707.3619v21.pdf>.
