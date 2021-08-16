# PriComp
A secure multi-party computation platform. This package has being tested under ubuntu 18.04 and debian 10.
## Compile PriComp
 1. Install cmake: ```apt install cmake```
 2. Install c++ compiler: ```apt install g++```
 3. Install gmp library: ```atp install libgmp-dev```
 4. Use cmake: ```cmake .```
    - To show debugging messages, set -DDEBUG=Y
    - To evaluate the performance. set -DEVAL_TIME=Y
 5. Run make: ```make```
 6. Two directories, bin and lib, will be created. in bin, you can try some secure protocols. sp_demo performs scalar product, stat_demo computes some statistics, arith_demo performs the four arithmetic operationsm and f_condition_demo demonstrate how to write secure conditional structure in Pricomp.
## Compile commodity server
After successfully making PriComp, you can change directory to commodity and run cmake and make to build the programs of the commodity server and the programs that can generate online (necessary for the commodity server) and offline random bits.
## Compile ML examples
After successfully making PriComp, you can change directory to ML and run cmake and make to build the programs which can conduct secure machine learning protocols.
    
