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
 * System_Model.h
 *
 *  Created on: Jul 5, 2016
 *      Author: Alexandra Nagy
 *              alex.nagy.phys@gmail.com
 *
 *      REFERENCES:
 *      - A. Nagy, V. Savona, Driven-dissipative quantum Monte Carlo method for open quantum systems, Physical Review A 97 (5), 052129, 2018
 *      - A. Nagy's PhD Dissertation, École polytechnique fédérale de Lausanne, Quantum Monte Carlo Approach to the non-equilibrium steady state of open quantum systems, 2020
 */

#include "Main.h"

class System_Model
{
public:

	string name;

	struct Model_Properties
		{
			int dimension, ySite, zSite, xSite;
			double J; //coupling constant
			double h; //magnetic field
			unsigned long int reference;
		};


	Model_Properties prop;

	int Get_Excitation(bitset<B> left, bitset<B> right);
	double DiagonalElement(bitset<B> left, bitset<B> right);
	double OffDiagonalElement(bitset<B> left, bitset<B> right);
	double Calc_Occ_List_Renorm(bitset<B> det, vector<int> & occ);


	double Analytical(int site);
	vector<vector<double>> MyMatrixProduct(vector<double> A, vector<double> b);
	vector<double> MyMatrixAtanSum(vector<vector<double>> a, vector<vector<double>> b, vector<double> Ii);

	void Init(string filename);
	void MPI_Data(MPI_Datatype & p, MPI_Datatype & s, MPI_Datatype & sp, MPI_Datatype & sh);
	void MPI_Free(MPI_Datatype & p, MPI_Datatype & s, MPI_Datatype & sp, MPI_Datatype & sh);
	Random_excitation Gen_Rnd_Excitation(gsl_rng * random_gen, bitset<B> & determinant, vector<int> occ, double & renorm);
	bitset<B> Connected_Orbs(bitset<B> det, int init, vector<int> & occ_c);
	double HamiltonianElement(bitset<B> left, bitset<B> right);
	string OutputGen(int cycle, int target_nbr, double shift, bitset<B> reference, double tau, string & paath, int seed, double shift_const);
};
