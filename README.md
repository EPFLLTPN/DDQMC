# DDQMC
This repository contains an optimized, highly parallel C++ implementation of the Driven Dissipative Quantum Monte Carlo (DDQMC) method for sampling the non-equilibrium steady state of many-body open quantum systems.
For more information, see:

- A. Nagy, V. Savona, Driven-dissipative quantum Monte Carlo method for open quantum systems, Physical Review A 97 (5), 052129, 2018
- A. Nagy's PhD Dissertation, École polytechnique fédérale de Lausanne, Quantum Monte Carlo Approach to the non-equilibrium steady state of open quantum systems, 2020

# Physical model
The present version calculates the magnetization for the 2-dimensional dissipative XYZ Heisenberg model in an external magnetic field. The code can easily be adjusted for different models and/or observables by redefining System_Model.h, System_Model.cpp

# Requirements

