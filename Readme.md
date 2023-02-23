# k Shortest Paths

A reimplementation of the code used for my research papers about parallel k shortest paths algorithms.

This includes:
- Data structures for
  - directed and undirected graphs
  - Shortest path trees
- Shortest path algorithms
  - Dijkstra (with early stopping)
  - DeltaStepping (without bucket hierarchy, but with early stopping)
- k shortest paths algorithms
  - Yens algorithm
  - Fengs algorithm
  - some mixed variants with additional heuristics described in my PhD thesis.
- Generators for testing
  - G(n,p)
  - Grid
- Tools
  - Parser for some file formats (Metis, edgelist)

This repository also includes data from the experiments described in my PhD thesis. The data can be found in `experiments/data` and needs to be decompressed (about 575 MB decompressed) in respective folders. Then the Jupyter Notebook files in `experiments/evaluation` can be used to regenerate the the plots from the thesis.

## Build

We use cmake and gcc 10.3 to build the code. For easy usage there is build script.

```
./build.sh <target> [<numKSPthreads>] [<numDSthreads>]
```

`<target>` is one of the following executables
- `AllTests`: runs all tests using gtest.
- `convert`: converts graph files.
- `generate`: generates graphs as metis files.
- `SingleThreadPerformance`: runs single thread performance experiments (mostly for quick testing).
- `PerformanceMeasurements`: runs performance experiments. Uses [<numKSPthreads>] and [<numDSthreads>] to run on a set number of threads.
- `CollectStatisticalData`: collects statistics about a k-shortest-path run. (path length, number of hops, size of shortest path sub trees, etc). This program does not measure runtime and runs only in a single thread to get the order of the data right.

The script creates a build folder named `build-release`. There all executables can be found.

## Usage

All executables have a CLI interface. Just run it the `-h` option or no option at all to get the help page printed out. Results of the experiments are printed to std::out while other information during runtime (if any) is printed to std::err. The output should be redirected into a file.

In order to run the tests, you first need to generate some graphs to run the tests on. The script `run-all-tests.sh` compiles the graph generator and tests, generates all graphs at the correct location (takes about 70 MB of disk space) and executes all tests.

## Open Source code used

- [CLI11](https://github.com/CLIUtils/CLI11) (single header version)
- Google Benchmark (git submodule)
- Google Test (git submodule)

## Cite

Please cite my PhD thesis when using this code in your reaserch project.

```
TODO Bibtex
```

## License

GNU GPL v3
