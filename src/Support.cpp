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
 * Support.cpp
 *
 *  Created on: Jul 5, 2016
 *      Author: Alexandra Nagy
 *              alex.nagy.phys@gmail.com
 *
 *      REFERENCES:
 *      - A. Nagy, V. Savona, Driven-dissipative quantum Monte Carlo method for open quantum systems, Physical Review A 97 (5), 052129, 2018
 *      - A. Nagy's PhD Dissertation, École polytechnique fédérale de Lausanne, Quantum Monte Carlo Approach to the non-equilibrium steady state of open quantum systems, 2020
 */

#include "Header/Main.h"

using namespace std;

int mapping (int nbr)
{
	return pow(nbr, 4);
}


string DivideComplex(double a, double b, double c, double d)
{
	stringstream ss;
	ss.str("");

	double real = (a*c + b*d)/(pow(c,2) + pow(d,2));
	double imag = (b*c - a*d)/(pow(c,2) + pow(d,2));

	if(imag > 0)
		ss << real <<"+"<<imag<<"i";
	else
		ss << real <<"-"<<abs(imag)<<"i";

	return ss.str();
}

vector<double> DivideComplexNbr(double a, double b, double c, double d)
{
	double real = (a*c + b*d)/(pow(c,2) + pow(d,2));
	double imag = (b*c - a*d)/(pow(c,2) + pow(d,2));

	vector<double> ret;

	ret.push_back(real);
	ret.push_back(imag);

	return ret;
}


int MPI_Hash(bitset<Bl> det, int range)
{
	int hash = 0;

	 boost::crc_optimal<32, 0x04C11DB7, 0xFFFFFFFF, 0xFFFFFFFF, true, true>  crc_ccitt2;

	 for(unsigned int index =0;index<det.size();++index)
	 {
	       crc_ccitt2(det[index]);
	 }

	 hash = crc_ccitt2();

	 hash = abs(hash%range);

	return hash;

}

void Hash_Test(int site, int range) {
	unsigned long int limit = pow(2, site);
	int val;
	string line;
	ofstream myfile;

	int count[range];
	for(int i=0; i<range;++i) count[i] =0;

	myfile.exceptions(ofstream::failbit | ofstream::badbit);
	try {
		myfile.open("/home/thean0/Desktop/CODE/src/test.txt");

		for (unsigned int i = 0; i < limit; ++i) {
cout << i << endl;
			bitset<Bl> state(i);

			val = MPI_Hash(state, range);
			count[val]++;
			myfile << i << "\t" << val << endl;
		}

		myfile.close();
	} catch (ofstream::failure e) {
		cout << "Exception opening/writing file";
	}

	for(int i=0; i<range;++i)
		cout << i << "\t" << count[i] << endl;
}


unsigned int split(const string &txt, vector<string> &strs, char ch)
{
	int pos = txt.find( ch );
	unsigned int initialPos = 0;
	strs.clear();

	// Decompose statement
	while( pos != (int)string::npos ) {
		strs.push_back( txt.substr( initialPos, pos - initialPos + 1 ) );
		initialPos = pos + 1;

		pos = txt.find( ch, initialPos );
	}

	// Add the last one
	strs.push_back(txt.substr(initialPos, min(pos, static_cast<int>(txt.size())) - initialPos ) );

	return strs.size();
}

