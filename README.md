# ERNN

## Purpose

This repository implements the algorithms mentioned in our VLDB 2024 paper (titled "Expanding Reverse Nearest Neighbors").

## Dataset

- All datasets used in the paper is downloaded from link http://www.diag.uniroma1.it/challenge9/download.shtml.
- The used datasets are available in the data folder of this repository (other oversized datasets can be obtained from link https://www.dropbox.com/scl/fo/9gmc1qjvbhiahjrk2ygsq/h?rlkey=fygkxf587t9dkgzcmw96efg2s&dl=0).
- The NH-d dataset is a small sample of the NH dataset that can be used to evaluate the exact algorithm.
- In addition, the test.tmp dataset was not used in the paper, but can be used for testing purposes.

## Compiling

- Use the following command to compile:
```
g++ -O3 -std=c++17 main.cpp -o run
```
- All algorithms were implemented in C++ and compiled using GNU GCC 8.5.0 with optimizations at the -O3 level. 
The experiments were conducted on a machine with an Intel Xeon 2.50 GHz CPU and 512 GB of memory, running 64-bit Red Hat Linux 8.5.0.

## Running

- Use the following command to run:
```
./run @1 @2 @3 @4 @5
```
  - Parameter @1 is the method name
    - Weight (select edges with the highest weights from the graph for weight reduction.)
    - Neighbor (select neighbors of the target vertex $f$ (if there are not enough neighbors of $f$, we choose neighbors of neighbors of $f$, and so on) and reduce the weights of the edges between vertex and the selected neighbors.)
    - Basic (this algorithm is described in Algorithm 2, which is a standard greedy algorithm.)
    - DBEI (this is our proposed algorithm, which is described in Algorithm 3. It uses a distance-based edge inspection technique to achieve early stopping using upper bounds.)
    - DEBIPLUS (DBEI combined with pruning strategy 1 (Lemma 4.12), which eliminates invalid edges.)
    - DEBISTAR (DBEIPLUS combined with pruning strategy 2 (Lemma 4.13), which eliminates dominated edges.)

  - Parameter @2 is the dataset name (stored in the data folder, without the .tmp extension)
  - Parameter @3 is the number of facility vertices, default is 1000
  - Parameter @4 is the number of rounds to run the experiment, where a target vertex is randomly selected for the experiment each time, default value is 50
  - Parameter @5 is the budget, i.e. the number of edges selected for the upgrade, default value is 4

## Running Example
- Execute the **DEBI** algorithm on the test.tmp dataset with 1000 vertices as facility vertices and perform 50 rounds of experiments with a budget of 4. 
```
./run DEBI test 1000 50 4
```
- Execute the **DEBISTAR** algorithm on the BAY.tmp dataset with 1000 vertices as facility vertices and perform 50 rounds of experiments with a budget of 4. 
```
./run DEBISTAR BAY 1000 50 4
```

## Interpretation of Results

- The following is the result of executing . /run DEBISTAR BAY 1000 50 4

> Round : 1, targetID = 152069

> Round : 2, targetID = 252656

> ...

> Round : 50, targetID = 153233

> Average Gain = 204.42, Average Time = 1.30738

- Round : 1, target = 152069 means that vertex 152069 is selected as the target facility to execute the DEBISTAR algorithm on BAY.tmp in round 1. (Due to randomness, the target vertex IDs selected during the actual run may be different.)
- Average Gain = 204.42 means that the average gain for 50 experiments is 204.42 and Average Time = 1.30738 means that the average run time is 1.30738 seconds.

## Additional Experiment
- Go to the index folder to verify Exp-9 in our paper (testing the use of indexing techniques).
  - (1) Compile by executing g++ -O3 -std=c++17 main.cpp -o run
  - (2) Run by executing ./run @1 @2 @3 @4
    - Parameter @1 is the dataset name (stored in the data folder, without the .tmp extension)
    - Parameter @2 is the number of facility vertices, default is 1000
    - Parameter @3 is the number of rounds to run the experiment, where a target vertex is randomly selected for the experiment each time, default value is 50
    - Parameter @4 is the random number seed (default value is 1)
- Example
```
./run test 1000 50 1
./run CT 1000 50 1
```
- No need to specify a budget as it is fixed at 1
- Note that this procedure will generate temporary indexes of significant size, so be aware of disk capacity when testing. Reset the fourth parameter (random number seed) for new rounds of randomness experiments when the disk capacity is insufficient.














