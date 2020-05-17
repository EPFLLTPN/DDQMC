# DDQMC
This repository contains an optimized, highly parallel C++ implementation of the Driven Dissipative Quantum Monte Carlo (DDQMC) method for sampling the non-equilibrium steady state of many-body open quantum systems.
For more information, see:

- A. Nagy, V. Savona, Driven-dissipative quantum Monte Carlo method for open quantum systems, Physical Review A 97 (5), 052129, 2018
- A. Nagy's PhD Dissertation, École polytechnique fédérale de Lausanne, Quantum Monte Carlo Approach to the non-equilibrium steady state of open quantum systems, 2020

## Physical model
The present version calculates the magnetization for the dissipative XYZ Heisenberg model in an external magnetic field (1D, 2D). The code can easily be adjusted for different models and/or observables by redefining System_Model.h, System_Model.cpp

## Some geography
Files are organized in the DDQMC repository as follows:
 - `src/` root directory of the program. It contains the main file (QMC_main.cpp), a support file (Support.cpp) that contains useful function definitions and a test unit (Test.cpp).
 - `src/Header/` contains the headers. 
 - `src/Classes/` contains the implementation files
 - `src/config.txt` the configuration file for the parameters of the simulation.
 - `src/model.txt` the model input file that defines the physical model and its parameters.
 - `result/` directory containing the result files.

## Requirements
### C++ compilers
DDQMC is written in C++. We have tested it using mpic++ compiler, that is a convenience wrappers for the underlying C++ compiler. Translation of an MPI program requires the linkage of the MPI-specific libraries which may not reside in one of the standard search directories of ld(1). It also often requires the inclusion of header files what may also not be found in a standard location. mpic++ passes its arguments to the underlying C++ compiler along with the -I, -L and -l options required by Open MPI programs.
### GSL - GNU Scientific Library
The GNU Scientific Library (GSL) is a numerical library for C and C++ programmers. It is free software under the GNU General Public License. GSL can be found in the gsl subdirectory on your nearest GNU mirror http://ftpmirror.gnu.org/gsl/.
### MPI
Since the parallelization is optimized for cluster computing, we used MVAPICH. The MVAPICH2 software, based on MPI 3.1 standard, delivers the best performance, scalability and fault tolerance for high-end computing systems and servers using InfiniBand, Omni-Path, Ethernet/iWARP, and RoCE networking technologies. 
### Boost
Boost provides free peer-reviewed portable C++ source libraries. Available from https://www.boost.org.

# Compilation
An example on how to compile the code is as follows:

Building file: ./Classes/QMC_state.cpp

`mpic++ -std=c++11 -I/gsl/2.1/include -I$BOOST_ROOT/include/ -O0 -g3 -Wall -c -o ./Classes/QMC_state.o ./Classes/QMC_state.cpp`

Building file: ./Classes/System_Model.cpp

`mpic++ -std=c++11 -I/gsl/2.1/include -I$BOOST_ROOT/include/ -O0 -g3 -Wall -c -o ./Classes/System_Model.o ./Classes/System_Model.cpp`

Building file: ./QMC_main.cpp

`mpic++ -std=c++11 -I/gsl/2.1/include -I$BOOST_ROOT/include/ -O0 -g3 -Wall -c -o ./QMC_main.o ./QMC_main.cpp`

Building file: ./Support.cpp

`mpic++ -std=c++11 -I/gsl/2.1/include -I$BOOST_ROOT/include/ -O0 -g3 -Wall -c -o ./Support.o ./Support.cpp`

Building target QMC

`mpic++ -o QMC ./Classes/QMC_state.o ./Classes/System_Model.o ./QMC_main.o ./Support.o -lgsl -L/gsl/2.1/lib -Xlinker -lgsl -lgslcblas -lm`

A more detailed explanation on how to compile parallel C++ code for a computational cluster can be found at https://scitas-data.epfl.ch/public/training/using_mpi.pdf.

# Usage
Example of a parallel run on a local machine:

`mpirun -np N QMC /path/to/config/ -`

where N is the number of parallel instances which needs to be EVEN. '/path/to/config' is the path to the folder where config.txt, model.txt and result/ is situated. The results of the simulations then will be printed in '/path/to/config/result/'. It is necessary that the folder exists.

Example of a parallel run file on a computational cluster using a SLURM workload manager:

`#!/bin/bash`

`#SBATCH --ntasks 80`

`#SBATCH --cpus-per-task 1`

`#SBATCH --time 10:00:00`

`#SBATCH --account=ltpn`

`module purge`

`module load gcc mvapich2 gsl boost`

`srun /path/to/QMC /path/to/config/ -`

# License
Except where otherwise stated, HANDE is under the 'Lesser GNU Public License v3.0'.
