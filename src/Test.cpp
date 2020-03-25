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
 * Test.cpp
 *
 *  Created on: Jul 5, 2016
 *      Author: Alexandra Nagy
 *              alex.nagy.phys@gmail.com
 *
 *      REFERENCES:
 *      - A. Nagy, V. Savona, Driven-dissipative quantum Monte Carlo method for open quantum systems, Physical Review A 97 (5), 052129, 2018
 *      - A. Nagy's PhD Dissertation, École polytechnique fédérale de Lausanne, Quantum Monte Carlo Approach to the non-equilibrium steady state of open quantum systems, 2020
 */

#include "Header/QMC_state.h"

void QMC_state::TestPerformance(Determinant d, int walk_nbr)
{
	int sp_re, sp_im;
	int ex;
	int proc_nbr;
	int n_sp;
	unsigned int * sp;
	double first, second;

	// Testing the for loop
	clock_t beg_uni = clock();
	for(int w=0; w<walk_nbr; ++w)
	{
		sp_re = 0;
		sp_im = 0;

		if(d.prop.conn_states.size() != 0)
		{
			ex = gsl_rng_uniform_int(gsl_random, d.prop.conn_states.size());

			sp_re = SpawningAttempt(d.prop.conn_ampl[ex], 1.0/d.prop.conn_states.size(), gsl_random, 1, 1, 1);
			sp_im = SpawningAttempt(d.prop.conn_ampl[ex], 1.0/d.prop.conn_states.size(), gsl_random, 1, 0, 1);

			if((sp_re != 0) || (sp_im != 0))
			{
				proc_nbr = MPI_Hash(d.prop.conn_states[ex], proc-1);
				spawned_send[proc_nbr].push_back({d.prop.conn_states[ex], sp_re, sp_im});
			}
		}
	}
	first =  (float( clock () - beg_uni ) /  CLOCKS_PER_SEC);
	//std::cout << "Execution time: \t" << (float( clock () - beg_uni ) /  CLOCKS_PER_SEC) << "  sec" << std::endl;
	for (int i = 0; i < proc - 1; ++i)
	{
		spawned_send[i].clear();
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Testing the for multinomial
	clock_t beg_multi = clock();

	if(d.prop.P_tot > 1)
		{
			n_sp = floor(d.prop.P_tot) + gsl_ran_binomial(gsl_random, d.prop.P_tot-floor(d.prop.P_tot), walk_nbr);
		}
		else
			n_sp = gsl_ran_binomial(gsl_random, d.prop.P_tot, walk_nbr);

		//MULTINOMIAL
		if(n_sp > 0)
		{
			sp = new unsigned int[d.prop.conn_prob.size()];
			gsl_ran_multinomial(gsl_random, d.prop.conn_prob.size(),n_sp, &d.prop.conn_prob[0], sp);

		for(unsigned int i = 0; i<d.prop.n_im;++i)
			if(sp[i] != 0){
				proc_nbr = MPI_Hash(d.prop.conn_states[i], proc - 1);
				spawned_send[proc_nbr].push_back({d.prop.conn_states[i], 0, ((d.prop.conn_ampl[i].imag())/(abs(d.prop.conn_ampl[i].imag())))*sp[i]});

			}

		for(unsigned int i = d.prop.n_im; i<d.prop.conn_states.size(); ++i)
			if(sp[i] != 0){
				proc_nbr = MPI_Hash(d.prop.conn_states[i], proc - 1);
				spawned_send[proc_nbr].push_back({d.prop.conn_states[i], ((d.prop.conn_ampl[i].imag())/(abs(d.prop.conn_ampl[i].imag())))*sp[i], 0});

		}
		}
		second = (float( clock () - beg_multi ) /  CLOCKS_PER_SEC) ;
		//std::cout << "Execution time: \t" << (float( clock () - beg_uni ) /  CLOCKS_PER_SEC) << "  sec" << std::endl;

		std::cout << walk_nbr <<"      "<<first<<"      "<<second<<"      "<<second-first<<std::endl;
}



