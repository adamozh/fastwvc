# FastWVC
A local search algorithm to solve the NP-Hard Minimum-Weighted Vertex Cover (MWVC) problem. This is my C++ implementation of the original FastWVC algorithm with some slight changes, referenced from the paper _Towards faster local search for minimum weight vertex cover on massive graphs_, published in 2018. A copy of the paper is uploaded in this repository. 

## Local Search
Instead of employing a typical deterministic universal, optimal or approximation algorithm, local search is a heuristic "stop anytime" method which moves between candidate solutions by applying local changes, until a candidate is deemed sufficiently optimal or until the time bound is reached. In this case, local changes are applied by deleting and removing candidate nodes from the current cover to reach a new candidate solution, until it has been running for 2 seconds. Details are commented in `fastwvc.cpp`.

## Running the algorithm
This algorithm was developed for and tested on the Kattis problem [Minimum Weighted Vertex Cover](https://open.kattis.com/problems/mwvc). It sets a time limit of 2 seconds. Specifics can be found here.
1. Build the executable with your favourite C++ compiler
```
clang++ fastwvc.cpp -std=c++20 -o fastwvc
```
2. Run the executable with the graph input
```
./fastwvc
8 9
1 1 999 1 1 1 999 100
0 1
1 2
1 4
2 3
2 5
3 6
4 5
5 6
6 7
```

## Results
This implementation scores a 35.69/40 on Kattis, where 40 means an optimal result for all 40 input graphs.
