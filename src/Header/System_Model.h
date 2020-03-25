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

	std::string name;

	std::unordered_map<std::bitset<Bl>, double>::const_iterator find_elm;

	struct Model_Properties
		{
			int xSite;
			double J_x, J_y, J_z, h, theta, gamma; //coupling constant
			double hop, double_hop;
			std::bitset<Bl> reference;
		};
	Model_Properties prop;

	double DiagonalElement(std::bitset<Bl> & left, int first);

	void Init(std::string filename);
	void MPI_Data(MPI_Datatype & p, MPI_Datatype & s, MPI_Datatype & sp, MPI_Datatype & sh, MPI_Datatype & bitset_type);
    void MPI_Free(MPI_Datatype & p, MPI_Datatype & s, MPI_Datatype & sp, MPI_Datatype & sh, MPI_Datatype & bitset_type);

	std::string OutputGen(int cycle, int target_nbr, double shift, std::bitset<Bl> reference, double tau, std::string & paath, int seed, double shift_const);
};
