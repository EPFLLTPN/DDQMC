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
 * Lindbladian_Test.h
 *
 *  Created on: Jul 5, 2016
 *      Author: Alexandra Nagy
 *              alex.nagy.phys@gmail.com
 *
 *      REFERENCES:
 *      - A. Nagy, V. Savona, Driven-dissipative quantum Monte Carlo method for open quantum systems, Physical Review A 97 (5), 052129, 2018
 *      - A. Nagy's PhD Dissertation, École polytechnique fédérale de Lausanne, Quantum Monte Carlo Approach to the non-equilibrium steady state of open quantum systems, 2020
 */

//For every different type of model, one has to rewrite System_Model.h and System_Model.cpp
#include "QMC_state.h"
#include <complex.h>
//#include "Main.h"

class Lindbladian_Test
{
public:
	std::vector<std::bitset<Bl>> states;
	int Nbr_of_States;
	int Nbr_of_Diagonal;
	double gamma;
	double time;

	double double_h, hop;

	struct Ampl
	{
		std::complex<double> a;
	};

	struct Ampl_Conn
	{
		std::complex<double> amplitude;
		std::bitset<Bl> st;
		std::vector<std::bitset<Bl>> row;
		std::vector<std::complex<double>> rowProb;
		std::vector<std::bitset<Bl>> col;
		std::vector<std::complex<double>> colProb;
		std::vector<std::bitset<Bl>> any;
		std::vector<std::complex<double>> anyProb;

		std::complex<double> diag;

		std::complex<double> Sz;

		int isdiagonal;
	};

	std::unordered_map<std::bitset<Bl>, Ampl_Conn> densityM;
	std::unordered_map<std::bitset<Bl>, Ampl> altered;
	std::unordered_map<std::bitset<Bl>, Ampl_Conn>::iterator find_elm;
	std::unordered_map<std::bitset<Bl>, Ampl_Conn>::iterator rho;
	std::unordered_map<std::bitset<Bl>, Ampl>::iterator alt;

	Lindbladian_Test();
	Lindbladian_Test(System_Model data, double tau);
	void Generate_Density_Matrix();
	void Generate_Connections(System_Model data, std::bitset<Bl> element);
	//int Spawn(double prob, double par_sgn);
	void ResetAltered();
	//int DeathCalculator(double jj, double ii, double gamma, int count, int parent, int child, int walker);
	void PrintDM();
};
