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
 * Statistics.cpp
 *
 *  Created on: Jul 6, 2016
 *      Author: Alexandra Nagy
 *              alex.nagy.phys@gmail.com
 *
 *      REFERENCES:
 *      - A. Nagy, V. Savona, Driven-dissipative quantum Monte Carlo method for open quantum systems, Physical Review A 97 (5), 052129, 2018
 *      - A. Nagy's PhD Dissertation, École polytechnique fédérale de Lausanne, Quantum Monte Carlo Approach to the non-equilibrium steady state of open quantum systems, 2020
 */

#include "../Header/Main.h"

void Statistics::Write_Balance_init(string filename) {

	string line;
	ofstream myfile;
	stringstream s;

	s << filename << "_BALANCE.time";

	myfile.exceptions(ofstream::failbit | ofstream::badbit);
	try {

		myfile.open(s.str().c_str());

		myfile << "Time step \t Master (" << master << ") \t ";

		for(int i=0; i<proc-1; ++i)
			myfile << "Slave (" << i << ") \t ";

		myfile << endl;
		myfile << "-----------------------------------------------------------" << endl;

		myfile.close();
	} catch (ofstream::failure e) {
		cout << "Exception opening/writing file";
	}
}

