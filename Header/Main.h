/* 
 * This file is part of the EPFLLTPN/DDQMC distribution (https://github.com/EPFLLTPN/DDQMC).
 * Copyright (c) 2016 Alexandra Nagy.
 * 
 * This program is free software: you can redistribute it and/or modify  
 * it under the terms of the GNU General Public License as published by  
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 * Main.h
 *
 *  Created on: Jul 5, 2016
 *      Author: Alexandra Nagy
 *              alex.nagy.phys@gmail.com
 *
 *      REFERENCES:
 *      - A. Nagy, V. Savona, Driven-dissipative quantum Monte Carlo method for open quantum systems, Physical Review A 97 (5), 052129, 2018
 *      - A. Nagy's PhD Dissertation, École polytechnique fédérale de Lausanne, Quantum Monte Carlo Approach to the non-equilibrium steady state of open quantum systems, 2020
 */

//Length of bitsets: CHECK IF YOUR SYSTEM IS SMALLER!!
#define B 150

#pragma once

#include <iostream>
#include <time.h>
#include <string>
#include <sstream>
#include <fstream>
#include <bitset>
#include <functional>
#include <stdlib.h>
#include <vector>
#include <unordered_map>
#include <bitset>
#include <math.h>
#include <cmath>
#include "mpi.h"
#include <gsl/gsl_rng.h>
#include <algorithm>


using namespace std;

//Global variables
extern string path;
extern int proc, myrank, master;
extern MPI_Request request;

//Helpful functions
int MPI_Hash(unsigned long int det, int range);
void Hash_Test(int site, int range);

//Structures
struct Pop_control
{
	int ctrl;
	int alert;
	int ref_on;
	double shift;
};

struct Spawned_data
{
	ulong det;
	int nbr;
};

struct Stats
{
	int N_w, N_w_old, det_nbr, E_div;
	double E;
};

struct Determinant_properties
{
	vector<int> occ_list;  //list of the occupied spin-orbitals (Heisenberg: spin-up sites)
	double renorm;   //renormalization for p_gen
};

struct Random_excitation
{
	double p_gen; // generation probability
	ulong child; //the generated child determinant
	double H_elm; //the Hamiltonian element between the 2 determinants (parent and child)
};

//Classes
class Statistics
{
public:

	void Write_Balance_init(string filename);
};





