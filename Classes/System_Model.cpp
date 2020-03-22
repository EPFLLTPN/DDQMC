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
 * System_model.cpp
 *
 *  Created on: Jul 5, 2016
 *      Author: Alexandra Nagy
 *              alex.nagy.phys@gmail.com
 *
 *      REFERENCES:
 *      - A. Nagy, V. Savona, Driven-dissipative quantum Monte Carlo method for open quantum systems, Physical Review A 97 (5), 052129, 2018
 *      - A. Nagy's PhD Dissertation, École polytechnique fédérale de Lausanne, Quantum Monte Carlo Approach to the non-equilibrium steady state of open quantum systems, 2020
 */

#include "../Header/System_Model.h"
#include <climits>


unsigned long toInt(std::string const &s);

void System_Model::Init(string filename)
{
	//basic properties initialization

	string line;
	ifstream myfile;

	myfile.exceptions ( ifstream::failbit | ifstream::badbit );
	try
	{
		prop.dimension = 1;
		prop.ySite = 0;
		prop.zSite = 0;

		myfile.open (filename.c_str());

		for(int i=0;i<3;++i) getline(myfile, line);
					name = line;

		name.erase(name.end()-1);

		for(int i=0;i<5;++i) getline(myfile, line);
			prop.xSite = atoi(line.c_str());

		for(int i=0;i<5;++i) getline(myfile, line);
			prop.J = atof(line.c_str());

		for(int i=0;i<5;++i) getline(myfile, line);
			prop.h = atof(line.c_str());

		for(int i=0;i<5;++i) getline(myfile, line);

		prop.reference = toInt(line);

		myfile.close();
	}
	catch (ifstream::failure e)
	{
		cout << "Exception opening/reading file";
	}
}

void System_Model::MPI_Data(MPI_Datatype & p, MPI_Datatype & s,  MPI_Datatype & sp, MPI_Datatype & sh)
{
	//propterties
	MPI_Datatype oldtypes[3];
	int blockcounts[3];
	MPI_Aint offsets[3], extent1, extent2;

	offsets[0] = 0;
	oldtypes[0] = MPI_INTEGER;
	blockcounts[0] = 4;

	MPI_Type_extent(MPI_INTEGER, &extent1);
	offsets[1] = 4 * extent1;
	oldtypes[1] = MPI_DOUBLE;
	blockcounts[1] = 2;

	MPI_Type_extent(MPI_DOUBLE, &extent2);
	offsets[2] = 4 * extent1 + 2*extent2;
	oldtypes[2] = MPI_UNSIGNED_LONG;
	blockcounts[2] = 1;

	MPI_Type_struct(3, blockcounts, offsets, oldtypes, &p);
	MPI_Type_commit(&p);


	///////////////////////////////////////////////////////////////////////////////////
	//statistics
	Stats proba;
	MPI_Datatype oldtypes2[2];
	int blockcounts2[2];
	MPI_Aint offsets2[2];

	blockcounts2[0] = 4;
	blockcounts2[1] = 1;

	oldtypes2[0] = MPI_INTEGER;
	oldtypes2[1] = MPI_DOUBLE;

	MPI_Address(&proba.N_w, &offsets2[0]);
	MPI_Address(&proba.E, &offsets2[1]);

	offsets2[1] = offsets2[1] - offsets2[0];
	offsets2[0] = 0;

	MPI_Type_struct(2, blockcounts2, offsets2, oldtypes2, &s);
	MPI_Type_commit(&s);

	///////////////////////////////////////////////////////////////////////////////////
	//spawned data
	Spawned_data sp_probe;
	MPI_Datatype oldtypes3[2];
	int blockcounts3[2];
	MPI_Aint offsets3[2];

	blockcounts3[0] = 1;
	blockcounts3[1] = 1;

	oldtypes3[0] = MPI_UNSIGNED_LONG;
	oldtypes3[1] = MPI_INTEGER;

	MPI_Address(&sp_probe.det, &offsets3[0]);
	MPI_Address(&sp_probe.nbr, &offsets3[1]);

	offsets3[1] = offsets3[1] - offsets3[0];
	offsets3[0] = 0;

	MPI_Type_struct(2, blockcounts3, offsets3, oldtypes3, &sp);
	MPI_Type_commit(&sp);


	///////////////////////////////////////////////////////////////////////////////////
	//population control
	Pop_control pop_probe;
	MPI_Datatype oldtypes4[2];
	int blockcounts4[2];
	MPI_Aint offsets4[2];

	blockcounts4[0] = 3;
	blockcounts4[1] = 1;

	oldtypes4[0] = MPI_INTEGER;
	oldtypes4[1] = MPI_DOUBLE;

	MPI_Address(&pop_probe.ctrl, &offsets4[0]);
	MPI_Address(&pop_probe.shift, &offsets4[1]);

	offsets4[1] = offsets4[1] - offsets4[0];
	offsets4[0] = 0;

	MPI_Type_struct(2, blockcounts4, offsets4, oldtypes4, &sh);
	MPI_Type_commit(&sh);

}

void System_Model::MPI_Free(MPI_Datatype & p, MPI_Datatype & s,  MPI_Datatype & sp, MPI_Datatype & sh)
{
	MPI_Type_free(&p);
	MPI_Type_free(&s);
	MPI_Type_free(&sp);
	MPI_Type_free(&sh);
}

Random_excitation System_Model::Gen_Rnd_Excitation(gsl_rng * random_gen, bitset<B> & determinant, vector<int> occ, double & renorm)
{
	//Choose a random connected excitation
	// do until find at least one spin-up with available excitations
	Random_excitation rnd_ex;
	vector<int> occ_child;
	bitset<B> virt_avail, child;
	rnd_ex.p_gen = 0;
	bool excit = false;
	int help, init, dest, n_spinup;
	n_spinup = occ.size();

	while((excit == false) || (occ.size() != 0))
	{
		//random selection of initial spin-up
		help = gsl_rng_uniform_int(random_gen, occ.size());
		init = occ[help];

		occ.erase(occ.begin() + help);

		//Does this have available excitation?
		//Connected_orbs gives the bitstring with 1-s at the available excitation sites
		//compl is the complement of the original determinant
		// AND operation gives the bitstring containing the available excitations
		occ_child.clear();
		virt_avail = ~determinant &= Connected_Orbs(determinant, init, occ_child);

		if(virt_avail.count() != 0)
		{
			excit = true;
			break;
		}
	}

	//if there is absolutely no available excitation from the determinant
	if(!excit)
		return rnd_ex;

	//Now we find the excitation which we choose from the available ones
	help = gsl_rng_uniform_int(random_gen, occ_child.size());
	dest = occ_child[help];

	//set the child determinant
	child = determinant;
	child[init] = 0;
	child[dest] = 1;
	rnd_ex.child = child.to_ulong();
	//set the matrix element
	rnd_ex.H_elm = - prop.J / 2.0;

	//Calculate the generation probability
	//P_gen=p(i)*p(f|i)*renorm
	//p(i)=1/nbr_spinup, p(f|i)=1/nbr_avail_from_init
	rnd_ex.p_gen = (1/(double)n_spinup)*(1/(double)occ_child.size())*renorm;

	return rnd_ex;
}

bitset<B> System_Model::Connected_Orbs(bitset<B> det, int init,  vector<int> & occ_c)
{
	if(init == (prop.xSite - 1))
	{
		if(det[init-1] == 0){
			det[init-1] = 1;
			occ_c.push_back(init-1);}
		else
			det[init-1] = 0;

		if(det[0] == 0){
			det[0] = 1;
			occ_c.push_back(0);}
		else
			det[0] = 0;
	}
	else if(init == 0)
	{
		if (det[init + 1] == 0){
			det[init + 1] = 1;
			occ_c.push_back(init+1);}
		else
			det[init + 1] = 0;

		if (det[(prop.xSite - 1)] == 0){
			det[(prop.xSite - 1)] = 1;
			occ_c.push_back((prop.xSite - 1));}
		else
			det[(prop.xSite - 1)] = 0;

	} else {
		if (det[init - 1] == 0){
			det[init - 1] = 1;
			occ_c.push_back(init-1);}
		else
			det[init - 1] = 0;

		if (det[init + 1] == 0){
			det[init + 1] = 1;
			occ_c.push_back(init+1);}
		else
			det[init + 1] = 0;
	}

	return det;
}

int System_Model::Get_Excitation(bitset<B> left, bitset<B> right)
{

	int counter = 0;

	for(int i=(prop.xSite-1); i>=0;--i)
	{
		counter += (left[i]^right[i]);
	}

	return counter;
}

double System_Model::DiagonalElement(bitset<B> left, bitset<B> right)
{
	/*Contribution to Hamiltonian from spin interactions:
	 * For a lattice the number of bonds stored as Nbonds.
	 * Bonds of type 0-0 or 1-1 will give a contribution of -J to the matrix element,
	 * 0-1, 1-0 bonds will give +J contribution.
	 *
	 * Contribution the Hamiltonian from external field:
	 * Each spin up gives a contribution of -h/2. Each spin down gives gives +h/2
	 */

		double element;
		int counter = 0;
		int Nbonds = 0;

		for(int i=(prop.xSite-1); i>=0;--i)
		{
			++Nbonds;

			if(i == (prop.xSite-1))
				counter += (left[0]^left[i]);
			else
				counter += (left[i+1]^left[i]);
		}

		element = - (prop.J / 4.0) * (Nbonds - 2.0 * counter);

		element = element - prop.h / 2.0 * (2.0 * left.count() - prop.xSite);

		return element;
}

double System_Model::OffDiagonalElement(bitset<B> left, bitset<B> right)
{
	double element;

	for(int i=(prop.xSite-1); i>0;--i)
	{
		if((left[i]^right[i]) == 1)
		{
			if(i == (prop.xSite -1))
			{
				if(left[0]^right[0] == 1 || left[i-1]^right[i-1])
				{
					element = - prop.J / 2.0;
					break;
				} else {
					element = 0;
					break;
				}
			}


			if((left[i-1]^right[i-1]) == 1)
			{
				element = - prop.J / 2.0;
				break;
			}
			else
			{
				element = 0;
				break;
			}
		}
	}

	return element;
}

double System_Model::HamiltonianElement(bitset<B> left, bitset<B> right)
{

	double element;

	//Check if the determinants are not different with more than one spin orbit
	//Calculate the corresponding matrix element
	switch(Get_Excitation(left, right))
	{
	case 0:
		element = DiagonalElement(left, right);
		break;
	case 2:
		element = OffDiagonalElement(left, right);
		break;
	default:
		element = 0;
		break;
	}

	return element;
}

double System_Model::Calc_Occ_List_Renorm(bitset<B> det, vector<int> & occ)
{
	int no_excit = 0;

	for (int i = (prop.xSite - 1); i >= 0; --i) {

		if (det[i] == 1) {
			occ.push_back(i);

			if(i == (prop.xSite - 1))
			{
				if (det[i - 1] == 1 && det[0] == 1)
				no_excit++;
			}
			else if(i == 0)
			{
				if (det[i + 1] == 1 && det[(prop.xSite - 1)]==1 )
				no_excit++;
			}
			else
			{
				if (det[i - 1] == 1 && det[i + 1] == 1)
				no_excit++;
			}
		}
	}

	return det.count()/(double)(det.count() - no_excit);
}

string System_Model::OutputGen(int cycle, int target_nbr, double shift, bitset<B> reference, double tau, string & paath, int seed, double shift_cst) {

	string line;
	ofstream myfile;
	stringstream s;
	time_t now = time(0);
	s << path << "result/" << name << "_J" << prop.J << "_site" << prop.xSite << "_Iter" << cycle << "_Nbr"<< target_nbr << "_S" << shift;
	paath = s.str();
	s << ".out";

	myfile.exceptions(ofstream::failbit | ofstream::badbit);
	try {

		myfile.open(s.str().c_str());

		myfile << "Model:\t" << name << endl;
		myfile << "J:\t" << prop.J << endl;
		myfile << "xSite:\t" << prop.xSite << endl;
		myfile << "Reference:\t";
			for(int i = (prop.xSite-1); i>=0; --i)
				myfile << reference[i];
		myfile << endl;
		myfile << "---------------------------" << endl;
		myfile << "Number of iteration:\t" << cycle << endl;
		myfile << "Target walker number:\t" << target_nbr << endl;
		myfile << "Initial energy shift:\t" << shift << endl;
		myfile << "Time step:\t" << tau << endl;
		myfile << "Random seed:\t" << seed << endl;
		myfile << "Population control damping:\t" << shift_cst << endl;
		myfile << "Analytical result\t" <<  (prop.xSite)/4.0 + (prop.xSite)* Analytical(prop.xSite) << endl;
		myfile << "Starting time:\t" << ctime(&now) << endl;
		myfile << "---------------------------" << endl;
		myfile << "Time \t Nbr of walk. \t Energy shift \t E \t Occ.det "<< endl;
		myfile << "--------------------------------------------" << endl;

		myfile.close();
	} catch (ofstream::failure e) {
		cout << "Exception opening/writing file";
	}

	return s.str();
}

vector<vector<double>> System_Model::MyMatrixProduct(vector<double> A, vector<double> b)
{
	vector<vector<double>> product;
	vector<double> row;

	for(int i = 0; i<A.size();++i)
	{
		for(int j=0; j<b.size();++j)
		{
			row.push_back(A[i]*b[j] / 2.0);
		}

		product.push_back(row);
		row.clear();
	}

	return product;
}

vector<double> System_Model::MyMatrixAtanSum(vector<vector<double>> a, vector<vector<double>> b, vector<double> Ii)
{
	vector<vector<double>> product;
	vector<double> row, result;
	double sum;

	for(int i = 0; i<a.size();++i)
		{
			for(int j=0; j<a[0].size();++j)
			{
				row.push_back(2.0 * atan(a[i][j]-b[i][j]));
			}

			product.push_back(row);
			row.clear();
		}

	for(int i =0; i<a[0].size(); ++i)
	{
		sum = 0;

		for(int j=0; j<a.size();++j)
		{
			sum += product[j][i];
		}

		sum /= (2.0 * (double)prop.xSite);
		sum += (Ii[i] * M_PI) / (double) prop.xSite;

		result.push_back(tan(sum));
	}

	return result;
}

double System_Model::Analytical(int site)
{
	int itNbr = 512;
	vector<double> I, zp, zz, ones;
	double E = 0;

	for(int i = 1; i < site; i +=2)
	{
		I.push_back(-(double)site/4.0 + i/2.0);
		zp.push_back(0.0);
		ones.push_back(1.0);
	}

	for(int i=0; i<itNbr; ++i)
	{
		zz = MyMatrixAtanSum(MyMatrixProduct(ones,zp),MyMatrixProduct(zp, ones), I);
		zp = zz;
	}

	for(int i = 0; i<zz.size(); ++i)
		E += (-2.0) / (1.0 + zz[i]*zz[i]);

	E = E*(-1.0)*prop.J/(double)site;

	return E;
}

unsigned long toInt(std::string const &s) {
    static const std::size_t MaxSize = CHAR_BIT*sizeof(unsigned long);
    if (s.size() > MaxSize) return 0; // handle error or just truncate?

    std::bitset<MaxSize> bits;
    std::istringstream is(s);
    is >> bits;
    return bits.to_ulong();
}
