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

int map (int nbr)
{
	return pow(nbr, 4);
}

void Hash_Collision(int site)
{
	unsigned long int limit = pow(2, site);
	vector<size_t> keys;
	int coll = 0;

		for (unsigned int i = 0; i < limit; ++i) {

			bitset<Bl> state(i);
			RepType st = BitsetToInt(state);

			 std::size_t seed = 0;
			 for (int i = 0; i <ByteS; ++i)
			     boost::hash_combine(seed, st[i] * 2654435761);

			 if(find(keys.begin(), keys.end(), seed) != keys.end())
				 ++coll;
			 else
				 keys.push_back(seed);

				cout << i <<"    "<<seed<< endl;

		}

		cout <<"COLLIEDED    :"<<coll << endl;
}


RepType BitsetToInt(std::bitset<Bl> state)
{
	RepType bytes;
	bitset<Bl> mask("1111111111111111");

	for(int i = 0; i<ByteS;++i)
	{
		bytes[i] = (unsigned int)((state >> (i*16)) & mask).to_ulong();
	}

	return bytes;
}

std::bitset<Bl> IntToBitset(RepType arr)
{
	std::bitset<Bl> state(0);

	for(int i = 0; i<ByteS;++i)
	{
		std::bitset<16> part(arr[i]);

		for(int j =0; j<16;++j)
			state[j + i*16] = part[j];
	}

	return state;
}



int IsDiagonal(std::bitset<Bl> & child, int & le)
{
	int count = 0;
	for(int N = 0; N<le;++N)
		if(child[N]^child[N+le])
				++count;

	if(count == 0)
		return 1;
	else
		return 0;
}

int MPI_Hash(RepType det, int range)
{
	int hash = 0;

	 boost::crc_optimal<32, 0x04C11DB7, 0xFFFFFFFF, 0xFFFFFFFF, true, true>  crc_ccitt2;

	 for(unsigned int index =0;index<ByteS;++index)
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
		myfile.open("/home/thean0/Desktop/OSW/src/test.txt");

		for (unsigned int i = 0; i < limit; ++i) {
cout << i << endl;
			bitset<Bl> state(i);
			RepType st = BitsetToInt(state);

			val = MPI_Hash(st, range);
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

