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

 If things doest not go wrong, two directories, bin and lib, will be created. In bin, you can try some secure protocols. sp_demo performs scalar product, stat_demo computes some statistics, arith_demo performs the four arithmetic operationsm and f_condition_demo demonstrate how to write secure conditional structure in Pricomp.
## How to run the examples
 Assume you have two servers S1 and S2 that have successfully compiled PriComp and want to play the role of Alice and Bob respectively. S1's IP address is 192.168.1.1, and S2's IP address is 192.168.1.2. 
 1. Go to directory conf and modify the cfg file in both S1 and S2
   - Run on different server: Set local_run=false
   - Use offline random bits: Set offline_commodity=true
   - Set offlineRNDDIR to the directory containining the files of offline random bits. The two files alice.dat and bob.dat in conf are the two files we generated. You can do "compile commodity server" and generate your own files.
   - Set IP in line 17. In this case, it should be: ```$host0, host1, host2 = (local_run) ? ["localhost", "localhost", "localhost"] : ["127.0.0.1", "192.168.1.1", "192.168.1.2"]```
   - Change the port numbers 30002 and 30003 to other port number if it's necessary.
   - Save the file (assume "example.cfg")
 2. Go to directory bin. You can find some execuable files. For example, arith_demo shows how to perform arithmetic operations. f_cond_demo shows how to securely compare two values. You can run these examples. Taking f_cond_demo as the example. If there are two parties S1 (Alice) and S2 (Bob). The two parties can execure the following command:
   - In S1: ```./f_cond_demo ../conf/example.cfg 1.0 1```
   - In S2: ```./f_cond_demo ../conf/example.cfg 3.0 2```

 To see the usage of these examples, run the executable files without giving any parameter.
## Compile commodity server
After successfully making PriComp, you can change directory to commodity and run cmake and make to build the programs of the commodity server and the programs that can generate online (necessary for the commodity server) and offline random bits.
## Compile ML examples
After successfully making PriComp, you can change directory to ML and run cmake and make to build the programs which can conduct secure machine learning protocols.
## Implement your own protocols
Before performing any secure operations PriComp provides, call SP_init first. When your protocol ends, call SP_clear. Examples can be found in the cpp files we provides in the root dircttory of priComp.

Currently PriComp only supports secure two-party computations. When you run prototols built using PriComp, you have to run Alice before Bob.
## Generate documents
You can generate PriComp API documents using doxygen. After you clone Pricomp in your computer:
  1. Install doxygen: ```apt install doxygen```
  2. Go to the root directory of PriComp and generate PriComp documents:```doxygen pricomp_doxygen.cfg```
  
A directory docs will be created. Open docs/index.html to read PriComp API documents.
