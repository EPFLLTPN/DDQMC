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

//Length of std::std::strings: CHECK IF YOUR SYSTEM IS SMALLER!!
#define Bl 100

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
#include <math.h>
#include <cmath>
#include "mpi.h"
#include <iomanip>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <algorithm>
#include <boost/crc.hpp>
#include <boost/functional/hash.hpp>
#include <boost/lexical_cast.hpp>

//Global variables
extern std::string path, restart_path;
extern int proc, myrank, master, rest_iter, iter, energy_freq;
extern MPI_Request request;
extern bool restart;
extern std::complex<double> im;

//Helpful functions
int MPI_Hash(std::bitset<Bl> det, int range);
void Hash_Test(int site, int range);
unsigned int split(const std::string &txt, std::vector<std::string> &strs, char ch);
std::string DivideComplex(double a, double b, double c, double d);
std::vector<double> DivideComplexNbr(double a, double b, double c, double d);

struct Determinant_properties
{
	std::complex<double> mx;
	std::complex<double> my;
	std::complex<double> mz;

	std::vector<std::bitset<Bl>> conn_states;
	std::vector<std::complex<double>> conn_ampl;
	std::vector<double> conn_prob;
	std::complex<double> diag_ampl;
	double P_tot;
	int n_im;
};

struct Determinant
{
	std::bitset<Bl> det;
	int diagonal;
	int re;
	int im;
	Determinant_properties prop;
};


struct Pop_control
{
	int ctrl;
	double shift;
};

struct Spawned_data
{
	std::bitset<Bl> det;
	int nbr_re;
	int nbr_im;
};

struct Stats
{
	//int d_r, d_i;
	int N_w, N_w_old, det_nbr, diag_walk_re, diag_walk_im, N_tot;
	std::complex<double> Mx, My, Mz, Div;
};




