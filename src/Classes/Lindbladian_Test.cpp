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
 * Lindbladian_Test.cpp
 *
 *  Created on: Jul 5, 2016
 *      Author: Alexandra Nagy
 *              alex.nagy.phys@gmail.com
 *
 *      REFERENCES:
 *      - A. Nagy, V. Savona, Driven-dissipative quantum Monte Carlo method for open quantum systems, Physical Review A 97 (5), 052129, 2018
 *      - A. Nagy's PhD Dissertation, École polytechnique fédérale de Lausanne, Quantum Monte Carlo Approach to the non-equilibrium steady state of open quantum systems, 2020
 */
#include "../Header/Lindbladian_Test.h"

using namespace std;

Lindbladian_Test::Lindbladian_Test(){}

Lindbladian_Test::Lindbladian_Test(System_Model data, double tau)
{
	//Calc properties
	Nbr_of_States = pow(2, data.prop.xSite) * pow(2, data.prop.xSite);
	Nbr_of_Diagonal = pow(2, data.prop.xSite);
	time = tau;
	gamma = data.prop.gamma;

	double_h = (data.prop.J_x-data.prop.J_y)/4.0;
	hop = (data.prop.J_x+data.prop.J_y)/4.0;

	//Create all the possible states, fill up the density matrix with zeros
	complex<double> c(0,0);
	Ampl a = {c};
	vector<complex<double>> hh;
	vector<bitset<Bl>> sg;
	vector<bitset<Bl>> em;
	int diag = 0;

	bitset<Bl> element;
	for(int i = 0; i<Nbr_of_Diagonal; ++i)
	{
		bitset<Bl> help(i);
		sg.push_back(help);
	}

	for(unsigned int i = 0; i<sg.size();++i)
	{
		for(unsigned int j = 0; j<sg.size(); ++j)
		{
			for(int first = 0 ; first < data.prop.xSite; ++first)
				element[data.prop.xSite - 1 - first] = sg[i][first];

			for(int second = 0 ; second < data.prop.xSite; ++second)
				element[(2*data.prop.xSite) - 1 - second] = sg[j][second];

			for(int n = (2* data.prop.xSite); n<Bl; ++n)
				element[n] = 0;

			if(sg[i] == sg[j])
				diag = 1;
			else
				diag = 0;

			Ampl_Conn ac = {c, element, em, hh, em, hh, em, hh, c, c, diag};
			densityM.insert({element, ac});
			altered.insert({element, a});
		}
	}


	//Calculate all the possible connections to every element
	unordered_map<std::bitset<Bl>, Ampl_Conn>::iterator e = densityM.begin();
	while(e != densityM.end())
	{
		Generate_Connections(data, e->second.st);


		cout << e->second.st;
		cout << "       Row" << endl;
		for(int i = 0; i<e->second.row.size();++i)
			cout <<"        "<<e->second.row[i] <<"        " << e->second.rowProb[i]<< endl;
		cout << "       Col" << endl;
		for(int i = 0; i<e->second.col.size();++i)
			cout <<"        "<<e->second.col[i] <<"        " << e->second.colProb[i]<< endl;
		cout << "       Any" << endl;
		for(int i = 0; i<e->second.any.size();++i)
			cout <<"        "<<e->second.any[i] <<"        " << e->second.anyProb[i]<< endl;
		cout << e->second.diag << endl;


		++e;
	}

	//put walkers on reference element
	find_elm = densityM.find(data.prop.reference);
	complex<double> cc(1,0);
	find_elm->second.amplitude = cc;

}

void Lindbladian_Test::PrintDM()
{
	complex<double> mz = 0;

	unordered_map<std::bitset<Bl>, Ampl_Conn>::iterator e = densityM.begin();
		while(e != densityM.end())
		{
			mz += e->second.Sz * e->second.amplitude;

			//if(e->second.isdiagonal)
				cout << e->second.st << "   "<< e->second.amplitude << endl;
			++e;
		}

	cout <<"mz:     " << mz <<endl;
}

void Lindbladian_Test::Generate_Density_Matrix()
{
	unordered_map<std::bitset<Bl>, Ampl_Conn>::iterator d = densityM.begin();
	while (d != densityM.end())
	{
		alt = altered.find(d->second.st);

		for(unsigned int r = 0; r<d->second.row.size();++r)
		{
			rho = densityM.find(d->second.row[r]);
			alt->second.a += d->second.rowProb[r] * rho->second.amplitude*time;
		}

		for(unsigned int c = 0; c<d->second.col.size();++c)
		{
			rho = densityM.find(d->second.col[c]);
			alt->second.a += d->second.colProb[c] * rho->second.amplitude*time;
		}

		for(unsigned int an = 0; an<d->second.any.size();++an)
		{
			rho = densityM.find(d->second.any[an]);
			alt->second.a +=  d->second.anyProb[an] * rho->second.amplitude*time;
		}

		alt->second.a += d->second.diag * d->second.amplitude*time;

		++d;
	}

	//Summing up the old with the altered
	unordered_map<std::bitset<Bl>, Ampl_Conn>::iterator rho_new = densityM.begin();
	while (rho_new != densityM.end())
	{
		alt = altered.find(rho_new->second.st);

		rho_new->second.amplitude += alt->second.a;

		++rho_new;
	}

	ResetAltered();
}

void Lindbladian_Test::Generate_Connections(System_Model data,  std::bitset<Bl> element)
{
	bitset<Bl> state_l;
	int left;
	int help_init;
	complex<double> im(0,1);

	find_elm = densityM.find(element);

	//COL
	for(int init = 0; init<data.prop.xSite;++init)
	{
		state_l = element;

		help_init = init;
		left = (data.prop.xSite + (help_init - 1) % data.prop.xSite) % data.prop.xSite;

		//////////////////////////////////////////////////////////////////////
		//Hopping type + Double hopping type
		state_l.flip(left);
		state_l.flip(init);

		find_elm->second.row.push_back(state_l);
		if(element[left] != element[init])
			find_elm->second.rowProb.push_back(im*hop);
		else
			find_elm->second.rowProb.push_back(im*double_h);
	}

	//ROW
	for(int init = data.prop.xSite; init<2*data.prop.xSite;++init)
	{
		state_l = element;

		help_init = init-data.prop.xSite;
		left = (data.prop.xSite + (help_init - 1) % data.prop.xSite) % data.prop.xSite;
		left += data.prop.xSite;

		//////////////////////////////////////////////////////////////////////
		//Hopping type + Double hopping type
		state_l.flip(left);
		state_l.flip(init);

		find_elm->second.col.push_back(state_l);
		if(element[left] != element[init])
			find_elm->second.colProb.push_back(-im*hop);
		else
			find_elm->second.colProb.push_back(-im*double_h);
	}

	//ANY
    for(int init = 0; init < data.prop.xSite; ++init)
    {
    	if((element[init] == 0) && (element[init + data.prop.xSite] == 0))
    	{
    	  	state_l = element;
    	  	state_l.flip(init);
    	  	state_l.flip(init + data.prop.xSite);

    	  	find_elm->second.any.push_back(state_l);
    	  	find_elm->second.anyProb.push_back(gamma);
    	}
    }

    //DIAGONAL
    double hii = data.DiagonalElement(find_elm->second.st, 0);
    double hjj = data.DiagonalElement(find_elm->second.st, 1);
    find_elm->second.diag = im*(hjj-hii)-0.5 * gamma * element.count();


    if(find_elm->second.isdiagonal)
    {
    find_elm->second.Sz = 0;
    for(int i = 0; i<data.prop.xSite; ++i)
    	find_elm->second.Sz += 0.5 * pow(-1.0, find_elm->second.st[i] + 1);

    find_elm->second.Sz /= (double)data.prop.xSite;
    }
    else
    	find_elm->second.Sz = 0;
}

/*
int Lindbladian_Test::Spawn(double prob, double par_sgn)
{
	int n_spawn = 0;
	double p_spawn = 0;
	double sign = 0;

	p_spawn = time * abs(prob);
	sign = - 1.0 * (prob)/abs(prob) * par_sgn;

	//Need to take into account that multiple spawn can appear
	n_spawn = floor(p_spawn);
	p_spawn -= n_spawn;

	if(p_spawn > gsl_rng_uniform(gsl_random))
		++n_spawn;

	//decide the sign of the spawned walkers
	if(n_spawn > 0)
		n_spawn = n_spawn * sign;

	return n_spawn;
}
*/
void Lindbladian_Test::ResetAltered()
{
	unordered_map<std::bitset<Bl>, Ampl>::iterator a = altered.begin();
	while(a != altered.end())
	{
		a->second.a = 0;
		++a;
	}
}

/*
int Lindbladian_Test::DeathCalculator(double jj, double ii, double gamma, int count, int parent, int child, int walker)
{
	int kill;
	double p_death, element_re, element_im;

	element_im = jj - ii;
	element_re = -0.5 * gamma * count;

	if(parent) //real parent
	{
		if(child)  //real child
			p_death = time * element_re * walker;
		else  //imaginary child
			p_death = time * element_im * walker;
	}
	else   //imaginary parent
	{
		if(child)  //real child
			p_death = time * element_im * walker;
		else  //imaginary child
			p_death = time * element_re * walker;
	}

	if(p_death > 0)
		kill = floor(p_death);
	else
		kill = ceil(p_death);

	p_death -= kill; //remaining chance

	if(abs(p_death) > gsl_rng_uniform(gsl_random))
	{
		if(p_death < 0)
			--kill;
		else
			++kill;
	}

	return kill;
}
*/
