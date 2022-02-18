# Umatch_paper
Data and code which generated test data and ran benchmarks for the paper "U-match factorization: sparse homological algebra, lazy cycle representatives, and dualities in persistent (co)homology." by Haibin Hang, Chad Giusti, Lori Ziegelmeier, and Gregory Henselman-Petrusek.

Preprint available at: https://arxiv.org/abs/2108.08831

To run the as described in the paper, you will need the Rust programming language installed. The following directions will work for most \*NIX-based operating system; they may require some alterations for other OSes.
1) Clone this repository
2) Navigate to the ```code``` subfolder of the repository and build the benchmarking code via: 
  ```cargo build --release```
3) Return to the base folder in the repositry.
4) To run a benchmark on a clique complex data set, use the following command:
  ```/usr/bin/time -v ./code/target/release/compare_decomp_clique ./clique/<DATA_SET_NAME> <DIMENSION> <FIELD> <PRINT_GENS> <BENCHMARK_TYPE>```
  
   The parameters are as follows:
   - `<DATA_SET_NAME>` is the name of one of the subfolders of ``./clique/``.
   - `<DIMENSION>` is the homological dimension in which to run the benchmark. 
    
     Note: If the folder does not contain a file called pairs_dim<DIMENSION>.csv, you will need to generate one using ```./code/target/release/save_clique_pairs```. 
  
   - `<FIELD>` is a prime number giving the order of the coefficient field
   - `<PRINT_GENS>` is a boolean (true/false) indicating whether to enumerate and print a list of homology generators.
   - `<BENCHMARK_TYPE>` takes values in the set {1, 2, 3, 4}.
  
     1 = use row operations and the complete boundary matrix
  
     2 = use row operations restricted to pivot pairs
  
     3 = use column operations and the complete boundary matrix
  
     4 = use column operations restricted to pivot pairs
    
    Benchmarks for 2D/3D cubical data use the same format, but require the use of ```./code/target/release/compare_decomp_cubical2d``` and ```./code/target/release/compare_decomp_cubcial3d```, respectively.
    

