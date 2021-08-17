# PriComp
A secure multi-party computation platform. PriComp has being tested under ubuntu 18.04 and debian 10. Pricomp should be able to run under any unix-like OS with cmake, c++ compiler, and gmp library.
## Compile PriComp
 1. Install cmake: ```apt install cmake```
 2. Install c++ compiler: ```apt install g++```
 3. Install gmp library: ```atp install libgmp-dev```
 4. Use cmake: ```cmake .```
    - To show debugging messages, set -DDEBUG=Y
    - To evaluate the performance. set -DEVAL_TIME=Y
 5. Run make: ```make```
 6. If things doest not go wrong, two directories, bin and lib, will be created. In bin, you can try some secure protocols. sp_demo performs scalar product, stat_demo computes some statistics, arith_demo performs the four arithmetic operationsm and f_condition_demo demonstrate how to write secure conditional structure in Pricomp.
## Compile commodity server
After successfully making PriComp, you can change directory to commodity and run cmake and make to build the programs of the commodity server and the programs that can generate online (necessary for the commodity server) and offline random bits.
## Compile ML examples
After successfully making PriComp, you can change directory to ML and run cmake and make to build the programs which can conduct secure machine learning protocols.
## Generate documents
You can generate PriComp API documents using doxygen:
  1. Install doxygen: ```apt install doxygen```
  2. Generate the document:```doxygen pricomp_doxygen.cfg```
