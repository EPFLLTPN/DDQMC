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
 * QMC_State.cpp
 *
 *  Created on: Jul 5, 2016
 *      Author: Alexandra Nagy
 *              alex.nagy.phys@gmail.com
 *
 *      REFERENCES:
 *      - A. Nagy, V. Savona, Driven-dissipative quantum Monte Carlo method for open quantum systems, Physical Review A 97 (5), 052129, 2018
 *      - A. Nagy's PhD Dissertation, École polytechnique fédérale de Lausanne, Quantum Monte Carlo Approach to the non-equilibrium steady state of open quantum systems, 2020
 */

#include "../Header/QMC_state.h"

using namespace std;

QMC_state::QMC_state() {
}

QMC_state::QMC_state(string filename, System_Model & sys) {
	string line;
	ifstream myfile;

	myfile.exceptions(ifstream::failbit | ifstream::badbit);
	try {
		myfile.open(filename.c_str());

		for (int i = 0; i < 3; ++i)	getline(myfile, line);
		this->pop_ctrl.shift = atof(line.c_str());

		for (int i = 0; i < 5; ++i)	getline(myfile, line);
		this->shift_cnst = atof(line.c_str());

		for (int i = 0; i < 5; ++i)	getline(myfile, line);
		this->tau = atof(line.c_str());

		for (int i = 0; i < 5; ++i)	getline(myfile, line);
		this->target_nbr = atoi(line.c_str());

		for (int i = 0; i < 5; ++i)	getline(myfile, line);
		this->cycle = atoi(line.c_str());

		for (int i = 0; i < 5; ++i)	getline(myfile, line);
		this->seed = atoi(line.c_str());

		for (int i = 0; i < 5; ++i)	getline(myfile, line);
		this->write_out_period = atoi(line.c_str());

		for (int i = 0; i < 5; ++i)	getline(myfile, line);
		energy_freq = atoi(line.c_str());

		for (int i = 0; i < 5; ++i)	getline(myfile, line);
		impS = exp(-atof(line.c_str()));

		for (int i = 0; i < 5; ++i)	getline(myfile, line);
		InitLimit = atoi(line.c_str());

		myfile.close();

	} catch (ifstream::failure e) {

		cout << "Exception opening/reading file" << endl;
	}

	gsl_rng_env_setup();
	gsl_random = gsl_rng_alloc(gsl_rng_default);
	gsl_rng_set(gsl_random, this->seed);

	//setting up the reference determinant and the starting one -> put it to the correct process
	this->reference = sys.prop.reference;

	ref_process = MPI_Hash(BitsetToInt(this->reference), proc-1);

	if(myrank == ref_process)
	{
		if(!restart)
			Generate_Determinant(BitsetToInt(this->reference), 150, 0, sys);
	}

	//initialize the temporary spawned walker storage vector
	spawned_send.resize(proc - 1);
	spawned_recv.resize(proc - 1);

	//variable initialization
	stats.N_diag = 0;
	stats.N_diag_old = 150;
	stats.N_tot = 0;
	stats.det_nbr = 0;
	stats.diag_det_nbr = 0;
	stats.diag_walk_re = 0;
	stats.diag_walk_im = 0;
	stats.Div = 0;
	stats.Mx = 0;
	stats.My = 0;
	stats.Mz = 0;

	//population control off
	pop_ctrl.ctrl = 0;

	ReqArray = new MPI_Request[proc-1];
			MPI_Request help;
			for(int i = 0; i<proc-1;++i)
				ReqArray[i] = help;
}


void QMC_state::Sl_DeterminantStep(System_Model & sys)
{
	int proc_nbr;
	int diag_re_re, diag_im_im;
	int diag_re_im, diag_im_re;
	int diag_count_re = 0;
	int diag_count_im = 0;
	int diag_det = 0;
	double parent_sign_re, parent_sign_im;
	int n_sp;
	unsigned int * sp;

	unordered_map<RepType, Determinant, ArrayHasher>::iterator d = determinant.begin();
	while(d != determinant.end())
	{

		// REAL PARENT
		if(abs(d->second.re) != 0)
		{
			parent_sign_re = d->second.re / abs(d->second.re);

			if (d->second.prop.P_tot > 1)
				n_sp = floor(d->second.prop.P_tot)*abs(d->second.re)+ gsl_ran_binomial(gsl_random,d->second.prop.P_tot - floor(d->second.prop.P_tot), abs(d->second.re));
			else
				n_sp = gsl_ran_binomial(gsl_random, d->second.prop.P_tot, abs(d->second.re));

			//MULTINOMIAL
			if (n_sp > 0)
			{
				sp = new unsigned int[d->second.prop.conn_prob.size()];
				gsl_ran_multinomial(gsl_random, d->second.prop.conn_prob.size(), n_sp,	&d->second.prop.conn_prob[0], sp);

				for(unsigned int i = 0; i<d->second.prop.conn_states.size(); ++i)
				{
					if (sp[i] != 0)
					{
						proc_nbr = MPI_Hash(d->second.prop.conn_states[i], proc - 1);

						if(abs(d->second.prop.conn_sign[i]) == 2)
						{
							if(d->second.initiator)
								spawned_send[proc_nbr].push_back({ d->second.prop.conn_states[i],0, parent_sign_re*(d->second.prop.conn_sign[i]/2)* sp[i], d->second.prop.conn_diag[i], 1});
							else
								spawned_send[proc_nbr].push_back({ d->second.prop.conn_states[i],0, parent_sign_re*(d->second.prop.conn_sign[i]/2)* sp[i], d->second.prop.conn_diag[i], 0});
						}
						else
						{
							if(d->second.initiator)
								spawned_send[proc_nbr].push_back({ d->second.prop.conn_states[i],parent_sign_re*d->second.prop.conn_sign[i]* sp[i], 0, d->second.prop.conn_diag[i], 1});
							else
								spawned_send[proc_nbr].push_back({ d->second.prop.conn_states[i],parent_sign_re*d->second.prop.conn_sign[i]* sp[i], 0, d->second.prop.conn_diag[i], 0});
						}
					}
				}

				delete[] sp;
			}
		}

		// IMAGINARY PARENT
		if(abs(d->second.im) != 0)
		{
			parent_sign_im = d->second.im / abs(d->second.im);

			if (d->second.prop.P_tot > 1)
				n_sp = floor(d->second.prop.P_tot)*abs(d->second.im)+ gsl_ran_binomial(gsl_random,d->second.prop.P_tot - floor(d->second.prop.P_tot), abs(d->second.im));
			else
				n_sp = gsl_ran_binomial(gsl_random, d->second.prop.P_tot, abs(d->second.im));

			//MULTINOMIAL
			if (n_sp > 0)
			{
				sp = new unsigned int[d->second.prop.conn_prob.size()];
				gsl_ran_multinomial(gsl_random, d->second.prop.conn_prob.size(), n_sp,	&d->second.prop.conn_prob[0], sp);

				for(unsigned int i = 0; i<d->second.prop.conn_states.size(); ++i)
				{
					if (sp[i] != 0)
					{
						proc_nbr = MPI_Hash(d->second.prop.conn_states[i], proc - 1);

						if(abs(d->second.prop.conn_sign[i]) == 2)
						{
							if(d->second.initiator)
								spawned_send[proc_nbr].push_back({ d->second.prop.conn_states[i],parent_sign_im*(-1.0)*(d->second.prop.conn_sign[i]/2)* sp[i],0, d->second.prop.conn_diag[i],1 });
							else
								spawned_send[proc_nbr].push_back({ d->second.prop.conn_states[i],parent_sign_im*(-1.0)*(d->second.prop.conn_sign[i]/2)* sp[i],0, d->second.prop.conn_diag[i],0 });
						}
						else
						{
							if(d->second.initiator)
								spawned_send[proc_nbr].push_back({ d->second.prop.conn_states[i],0,parent_sign_im*d->second.prop.conn_sign[i]* sp[i], d->second.prop.conn_diag[i],1 });
							else
								spawned_send[proc_nbr].push_back({ d->second.prop.conn_states[i],0,parent_sign_im*d->second.prop.conn_sign[i]* sp[i], d->second.prop.conn_diag[i],0 });
						}
					}
				}

				delete[] sp;
			}
		}

		//Cloning/killing step
		//dying probability times the whole population of the determinant
		//number that dies for sure plus stochastic part
		complex<double> d_a;
		if(d->second.diagonal)
		{
			//d_a = d->second.prop.diag_ampl - pop_ctrl.shift;
			stats.N_diag += abs(d->second.re) + abs(d->second.im);
		}
		//else
			//d_a = d->second.prop.diag_ampl;
		d_a = d->second.prop.diag_ampl - pop_ctrl.shift;

		diag_re_im = Death(d_a, 1, 0, abs(d->second.re), parent_sign_re);
		diag_im_re = Death(d_a, 0, 1, abs(d->second.im), parent_sign_im);
		diag_re_re = Death(d_a, 1, 1, abs(d->second.re), parent_sign_re);
		diag_im_im = Death(d_a, 0, 0, abs(d->second.im), parent_sign_im);


		d->second.re += diag_re_re;
		d->second.re += diag_im_re;
		d->second.im += diag_im_im;
		d->second.im += diag_re_im;


		stats.N_tot += abs(d->second.re) + abs(d->second.im);

		if((abs(d->second.re) != 0) || (abs(d->second.im) != 0))
		{
			if (iter % energy_freq == 0)
			{
				complex<double> rho((double)d->second.re, (double)d->second.im);

				stats.Mx += rho * d->second.prop.mx;
				stats.My += rho * d->second.prop.my;
				stats.Mz += rho * d->second.prop.mz;

				if(d->second.diagonal)
					stats.Div += rho;
			}

			if(d->second.diagonal){
				diag_count_re += abs(d->second.re);
				diag_count_im += abs(d->second.im);
				diag_det += 1;
			}
			else
			{
				if((abs(d->second.re)+abs(d->second.im)) >= InitLimit)
					d->second.initiator = 1;
				else
					d->second.initiator = 0;
			}

			++d;
		}
		else
		{
			d = determinant.erase(d);
		}
	}

	stats.diag_walk_re = diag_count_re;
	stats.diag_walk_im = diag_count_im;
	stats.det_nbr = determinant.size();
	stats.diag_det_nbr = diag_det;
}

void QMC_state::Sl_SendStat(MPI_Datatype & stat_type)
{

	MPI_Isend(&stats, 1, stat_type, master, 0, MPI_COMM_WORLD, &request);
}

void QMC_state::PrintDensity(complex<double> trace)
{
	unordered_map<RepType, Determinant, ArrayHasher>::iterator d = determinant.begin();
	while(d != determinant.end())
	{
		if(d->second.diagonal == 1)
		{
			complex<double> elm(d->second.re, d->second.im);
			cout << IntToBitset(d->second.det) <<"   "<< (elm/abs(trace)).real()  << "    " << (elm/abs(trace)).imag() << endl;
		}
		++d;
	}
}

void QMC_state::Ms_RecvStat(MPI_Datatype & stat_type)
{
	Stats stat_temp;

	for(int i = 0; i < proc - 1; ++i)
	{
		MPI_Recv(&stat_temp, 1, stat_type, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		stats.N_diag += stat_temp.N_diag;
		stats.N_tot += stat_temp.N_tot;
		stats.det_nbr += stat_temp.det_nbr;
		stats.diag_det_nbr += stat_temp.diag_det_nbr;
		stats.diag_walk_re += stat_temp.diag_walk_re;
		stats.diag_walk_im += stat_temp.diag_walk_im;
		stats.Mx += stat_temp.Mx;
		stats.My += stat_temp.My;
		stats.Mz += stat_temp.Mz;
		stats.Div += stat_temp.Div;
	}

}


int QMC_state::Death(std::complex<double> & diag_ampl, int parent, int child, int walker, int par_sign)
{
	int kill;
	double p_death;
	int sign;

	//Cloning/killing step
	//dying probability times the whole population of the determinant
	//number that dies for sure plus stochastic part

	//Calculate the probability that spawning is successful
	if(parent) //real parent
	{
		if(child){ //trying to spawn real child
			p_death = tau * abs(diag_ampl.real());
			sign =  (diag_ampl.real())/abs(diag_ampl.real()) * par_sign;}
		else{   //trying to spawn imaginary child
			p_death = tau * abs(diag_ampl.imag());
			sign =  (diag_ampl.imag())/abs(diag_ampl.imag()) * par_sign;}
	}
	else   //imaginary parent
	{
		if(child){ //trying to spawn real child
			p_death = tau * abs(diag_ampl.imag());
			sign =  (-1.0)*(diag_ampl.imag())/abs(diag_ampl.imag()) * par_sign;}
		else{   //trying to spawn imaginary child
			p_death = tau * abs(diag_ampl.real());
			sign =  (diag_ampl.real())/abs(diag_ampl.real()) * par_sign;}
	}


	if(p_death > 1)
	{
		kill = floor(p_death)*walker+ gsl_ran_binomial(gsl_random,p_death - floor(p_death),walker);
	} else
		kill = gsl_ran_binomial(gsl_random, p_death,walker);

	if (kill > 0)
		kill = kill * sign;

	return kill;
}

void QMC_state::Ms_CalcShift()
{
	//shift calculation
		pop_ctrl.shift = pop_ctrl.shift + (shift_cnst / tau) * log10((double)stats.N_diag / (double)stats.N_diag_old);

}

void QMC_state::Ms_SendShift(MPI_Datatype & pop_type)
{

	for(int i = 0; i < proc - 1; ++i)
	{
		MPI_Isend(&pop_ctrl, 1, pop_type, i, 4, MPI_COMM_WORLD, &request);
	}

}

void QMC_state::Sl_RecvShift(MPI_Datatype & pop_type)
{
	MPI_Recv(&pop_ctrl, 1, pop_type, master, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

void QMC_state::Sl_SendSpawned(MPI_Datatype & sp_type)
{
	for(int i=0; i<proc - 1; ++i)
		MPI_Isend(&spawned_send[i][0], spawned_send[i].size(), sp_type, i, 2, MPI_COMM_WORLD, &ReqArray[i]);
}


void QMC_state::Sl_RecvSpawned(MPI_Datatype & sp_type,  MPI_Status & status, System_Model & sys)
{
	int length;

	for(int i = 0; i < proc - 1; ++i)
	{
		MPI_Probe(MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, &status);
		MPI_Get_count(&status, sp_type, &length);
		spawned_recv[status.MPI_SOURCE].resize(length);
		MPI_Recv(&spawned_recv[status.MPI_SOURCE][0], length, sp_type, status.MPI_SOURCE, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		Sl_SortSpawned(status.MPI_SOURCE, sys);
	}
}


void QMC_state::Sl_SortSpawned(int proc, System_Model & sys)
{
	for(unsigned int i = 0; i < spawned_recv[proc].size(); ++i)
	{

		find_det = determinant.find(spawned_recv[proc][i].det);

		if(find_det == determinant.end())
		{
			if(spawned_recv[proc][i].init_flag || spawned_recv[proc][i].diag)
				Generate_Determinant(spawned_recv[proc][i].det, spawned_recv[proc][i].nbr_re, spawned_recv[proc][i].nbr_im, sys);
		}
		else
		{
			find_det->second.re += spawned_recv[proc][i].nbr_re;
			find_det->second.im += spawned_recv[proc][i].nbr_im;
		}
	}
}


void QMC_state::Generate_Determinant(RepType determ, int nbr_re, int nbr_im, System_Model & sys)
{
	int count;
	vector<int> where;
	Determinant d;

	d.det = determ;
	d.re = nbr_re;
	d.im = nbr_im;

	if((abs(nbr_re) + abs(nbr_im)) >= InitLimit)
		d.initiator = 1;
	else
		d.initiator = 0;

	Determinant_properties properties;
	properties.P_tot = 0;

	bitset<Bl> det = IntToBitset(determ);

	////////////////////////////////////////////////////////////////////////////////
	//Calculating estimator element for the given density matrix state
	count = 0;
	for(int N = 0; N<sys.prop.Nsite;++N)
		if(det[N]^det[N+sys.prop.Nsite]){
			++count;
			where.push_back(N);}

	//DIAGONAL
	if(count == 0)
	{
		d.diagonal = 1;
		for(int i = 0; i<sys.prop.Nsite; ++i)
			properties.mz += pow(-1.0, det[i] + 1);

		properties.mz /= (double)sys.prop.Nsite;
	}
	//OFF-DIAGONAL
	else if(count == 1)
	{
		d.diagonal = 0;

		properties.mx = 1.0;

		if(det[where[0]] == 1)
			properties.my = -im;
		else
			properties.my = im;

		properties.mx /= (double)sys.prop.Nsite;
		properties.my /= (double)sys.prop.Nsite;
	}
	else
		d.diagonal = 0;

	////////////////////////////////////////////////////////////////////////////////
	//Importance sampling weight factors
	if(d.diagonal)
	{
		d.weight = 1.0;
		d.initiator = 1;
	}
	else
		d.weight = impS;


	properties.mz *= 1/d.weight;
	properties.my *= 1/d.weight;
	properties.mx *= 1/d.weight;

	////////////////////////////////////////////////////////////////////////////////
	//Filling up all the possible connections
	bitset<Bl> state_r, state_d, state_field;
	double child_weight;
	complex<double> hh, dhh, aa, field;
	int pos;
	int col_right, row_down;
	int rig, down;

	hh = im*sys.prop.hop*(1/d.weight);
	dhh = im*sys.prop.double_hop*(1/d.weight);
	aa = sys.prop.gamma*(1/d.weight);
	field = im*sys.prop.h*(1/d.weight);

	//COL
	for(int row = 0; row < sys.prop.ySite; ++row)
	{
		row_down = (sys.prop.ySite + (row+1)%sys.prop.ySite)%sys.prop.ySite;

		for(int col = 0; col < sys.prop.xSite; ++col)
		{
			//PREPARATION
			col_right = (sys.prop.xSite + (col+1)%sys.prop.xSite)%sys.prop.xSite;
			rig = row*sys.prop.xSite + col_right;
			down =  row_down*sys.prop.xSite + col;
			pos = row*sys.prop.xSite + col;

			state_r = det;
			state_d = det;

			if(sys.prop.h != 0)
			{
				state_field = det;
				Generate_Field_Col(properties, state_field, pos, sys, field, 1);
			}

			//////////////////////////////////////////////////////////////////////
			//Hopping type + Double hopping type
			Generate_Connection_Col(properties, det, state_r, pos, rig, sys, hh, dhh);
			Generate_Connection_Col(properties, det, state_d, pos, down, sys, hh, dhh);

		}
	}

	//ROW
	for(int row = 0; row < sys.prop.ySite; ++row)
	{
		row_down = (sys.prop.ySite + (row+1)%sys.prop.ySite)%sys.prop.ySite;

		for(int col = 0; col < sys.prop.xSite; ++col)
		{
			//PREPARATION
			col_right = (sys.prop.xSite + (col+1)%sys.prop.xSite)%sys.prop.xSite;
			rig = row*sys.prop.xSite + col_right + sys.prop.Nsite;
			down =  row_down*sys.prop.xSite + col + sys.prop.Nsite;
			pos = row*sys.prop.xSite + col + sys.prop.Nsite;

			state_r = det;
			state_d = det;

			if(sys.prop.h != 0)
			{
				state_field = det;
				Generate_Field_Col(properties, state_field, pos, sys, field, 0);
			}

			//////////////////////////////////////////////////////////////////////
			//Hopping type + Double hopping type
			Generate_Connection_Row(properties, det, state_r, pos, rig, sys, hh, dhh);
			Generate_Connection_Row(properties, det, state_d, pos, down, sys, hh, dhh);
		}
	}

	//ANY
    for(int init = 0; init < sys.prop.Nsite; ++init)
    {
    	if((det[init] == 1) && (det[init + sys.prop.Nsite] == 1))
    	{
    	  	state_r = det;
    	  	state_r.flip(init);
    	  	state_r.flip(init + sys.prop.Nsite);
    		if(IsDiagonal(state_r, sys.prop.Nsite))
    			child_weight = 1.0;
    		else
    			child_weight = impS;

    	  	properties.conn_states.push_back(BitsetToInt(state_r));
    	  	properties.conn_sign.push_back(int((aa.real()*child_weight)/abs(aa.real()*child_weight)));
    	  	properties.conn_prob.push_back(abs(aa*child_weight)*tau);
    	  	properties.P_tot += abs(aa*child_weight)*tau;

			if(child_weight == 1.0)
				properties.conn_diag.push_back(1);
			else
				properties.conn_diag.push_back(0);
    	}
    }

    for(unsigned int i = 0; i<properties.conn_prob.size(); ++i)
    	properties.conn_prob[i] /= properties.P_tot;


    double h_ii = sys.DiagonalElement(det, 0);
    double h_jj = sys.DiagonalElement(det, 1);
    properties.diag_ampl = im*(h_jj-h_ii)-0.5 * sys.prop.gamma * det.count();

	d.prop = properties;


//	if(properties.conn_prob.size() != properties.n_im)
	//{
/*
	bitset<Bl> alll = IntToBitset(d.det);
	bitset<9> iii;
	bitset<9> jjj;

	for(int t = 0; t<9;++t)
	{
		iii[t] = alll[t + 9];
		jjj[t] = alll[t];
	}


	cout << iii.to_ullong() << "  " << jjj.to_ullong() <<endl;
	cout << "       All" << endl;
	for(unsigned int i = 0; i<d.prop.conn_states.size();++i)
	{
		bitset<Bl> all = IntToBitset(d.prop.conn_states[i]);
		bitset<9> ii;
		bitset<9> jj;

			for(int t = 0; t<9;++t)
			{
				ii[t] = all[t + 9];
				jj[t] = all[t];
			}

	//	cout <<"        "<< ii.to_ullong() << "  " << jj.to_ullong() <<"        " << d.prop.conn_prob[i] <<"        " << d.prop.conn_prob[i]<<"         "<< d.prop.diag_ampl<< endl;
		cout <<"        "<< IntToBitset(d.prop.conn_states[i]) <<"        " << d.prop.conn_prob[i] <<"        " << d.prop.conn_sign[i]<<"         "<< d.prop.diag_ampl<< endl;

	}*/
//	}

	determinant.insert({determ, d});
}


void QMC_state::Generate_Field_Col(Determinant_properties &properties, std::bitset<Bl> &st, int & pos, System_Model & sys, std::complex<double> &ff, int col)
{
	double child_w;

	//////////////////////////////////////////////////////////////////////
	//Hopping type + Double hopping type
	st.flip(pos);

	if (IsDiagonal(st, sys.prop.Nsite))
		child_w = 1.0;
	else
		child_w = impS;

	properties.conn_states.push_back(BitsetToInt(st));

		//FIELD IN DIRECTION X
	if (sys.prop.theta == 0)
	{
		if (col) {
			properties.conn_sign.push_back(2*int((child_w * ff.imag()) / abs(child_w * ff.imag())));
			properties.conn_prob.push_back(abs(ff * child_w) * tau);
			properties.P_tot += abs(ff * child_w) * tau;
		} else {
			properties.conn_sign.push_back(2*int((-child_w * ff.imag()) / abs(-child_w * ff.imag())));
			properties.conn_prob.push_back(abs(-ff * child_w) * tau);
			properties.P_tot += abs(-ff * child_w) * tau;
		}
	}
	else //FIELD IN Y DIRECTION
	{
		if(col)
		{
			if(st[pos] == 0) //sigma+
			{
				properties.conn_sign.push_back(int((child_w * ff.imag()) / abs(child_w * ff.imag())));
				properties.conn_prob.push_back(abs(ff * child_w) * tau);
				properties.P_tot += abs(ff * child_w) * tau;
			}
			else //sigma-
			{
				properties.conn_sign.push_back(int((-child_w * ff.imag()) / abs(-child_w * ff.imag())));
				properties.conn_prob.push_back(abs(-ff * child_w) * tau);
				properties.P_tot += abs(-ff * child_w) * tau;
			}
		}
		else
		{
			if(st[pos] == 1) //sigma+
			{
				properties.conn_sign.push_back(int((-child_w * ff.imag()) / abs(child_w * ff.imag())));
				properties.conn_prob.push_back(abs(-ff * child_w) * tau);
				properties.P_tot += abs(-ff * child_w) * tau;
			}
			else //sigma-
			{
				properties.conn_sign.push_back(int((child_w * ff.imag()) / abs(child_w * ff.imag())));
				properties.conn_prob.push_back(abs(ff * child_w) * tau);
				properties.P_tot += abs(ff * child_w) * tau;
			}
		}
	}

	if (child_w == 1.0)
		properties.conn_diag.push_back(1);
	else
		properties.conn_diag.push_back(0);
}

void QMC_state::Generate_Connection_Col(Determinant_properties &properties,bitset<Bl> &det, bitset<Bl> &st, int &pos, int &ch_pos, System_Model & sys,complex<double> &hh, complex<double> &dhh)
{
	double child_w;

 	//////////////////////////////////////////////////////////////////////
	//Hopping type + Double hopping type
	st.flip(pos);
	st.flip(ch_pos);

	if (IsDiagonal(st, sys.prop.Nsite))
		child_w = 1.0;
	else
		child_w = impS;

	if(det[pos] != det[ch_pos]){
		properties.conn_states.push_back(BitsetToInt(st));

		properties.conn_sign.push_back(2*int((child_w*hh.imag())/abs(child_w*hh.imag())));
		properties.conn_prob.push_back(abs(hh*child_w)*tau);
		properties.P_tot += abs(hh*child_w)*tau;

		if(child_w == 1.0)
			properties.conn_diag.push_back(1);
		else
			properties.conn_diag.push_back(0);
	}
	else{
		if(sys.prop.double_hop != 0)
		{
			properties.conn_states.push_back(BitsetToInt(st));

			properties.conn_sign.push_back(2*int((child_w*dhh.imag())/abs(child_w*dhh.imag())));
			properties.conn_prob.push_back(abs(dhh*child_w)*tau);
			properties.P_tot += abs(dhh*child_w)*tau;

			if(child_w == 1.0)
				properties.conn_diag.push_back(1);
			else
				properties.conn_diag.push_back(0);
		}
	}
}

void QMC_state::Generate_Connection_Row(Determinant_properties &properties,bitset<Bl> &det, bitset<Bl> &st, int &pos, int &ch_pos, System_Model & sys,complex<double> &hh, complex<double> &dhh)
{
	double child_w;

	//////////////////////////////////////////////////////////////////////
	//Hopping type + Double hopping type
	st.flip(pos);
	st.flip(ch_pos);

	if (IsDiagonal(st, sys.prop.Nsite))
		child_w = 1.0;
	else
		child_w = impS;

	if(det[pos] != det[ch_pos]){
		properties.conn_states.push_back(BitsetToInt(st));

		properties.conn_sign.push_back(2*int((-child_w*hh.imag())/abs(-child_w*hh.imag())));
		properties.conn_prob.push_back(abs(-hh*child_w)*tau);
		properties.P_tot += abs(-hh*child_w)*tau;

		if(child_w == 1.0)
			properties.conn_diag.push_back(1);
		else
			properties.conn_diag.push_back(0);
	}
	else{
		if(sys.prop.double_hop != 0)
		{
			properties.conn_states.push_back(BitsetToInt(st));

			properties.conn_sign.push_back(2*int((-child_w*dhh.imag())/abs(-child_w*dhh.imag())));
			properties.conn_prob.push_back(abs(-dhh*child_w)*tau);
			properties.P_tot += abs(-dhh*child_w)*tau;

			if(child_w == 1.0)
				properties.conn_diag.push_back(1);
			else
				properties.conn_diag.push_back(0);
		}
	}
}


void QMC_state::Cleaner()
{
	stats.N_diag_old = stats.N_diag;
	stats.N_diag = 0;
	stats.N_tot = 0;
	stats.det_nbr = 0;
	stats.diag_det_nbr = 0;
	stats.diag_walk_re = 0;
	stats.diag_walk_im = 0;
	stats.Div = 0;
	stats.Mx = 0;
	stats.My = 0;
	stats.Mz = 0;

	for (int i = 0; i < proc - 1; ++i)
	{
		vector<Spawned_data>().swap(spawned_send[i]);
		vector<Spawned_data>().swap(spawned_recv[i]);
	}

	MPI_Request help;
	for(int i = 0; i<proc-1;++i)
		ReqArray[i] = help;
}

void QMC_state::Ms_WriteStat(string filename, int time)
{
	ofstream myfile;
	stringstream s;

	myfile.exceptions(ofstream::failbit | ofstream::badbit);
	try {

		myfile.open(filename.c_str(), ios::app);

		if(iter%energy_freq == 0)
		{
			myfile << time <<" \t " << stats.N_tot <<" \t " << stats.N_diag << "\t" <<
			pop_ctrl.shift <<" \t " << stats.det_nbr << "\t" <<
			stats.Mx / stats.Div  << "\t" <<
			stats.My / stats.Div  << "\t" <<
			stats.Mz / stats.Div << "\t" <<
			stats.diag_walk_re << "\t" << stats.diag_walk_im << "\t" <<
			(stats.Mz/stats.Div).real() << "\t" << (stats.Mz/stats.Div).imag() << endl;
		}
		else
			myfile << time <<" \t " << stats.N_diag << "\t"<< stats.diag_walk_re + stats.diag_walk_im << "\t" << pop_ctrl.shift <<" \t " << stats.det_nbr << "\t NaN \t NaN \t NaN \t"<< stats.diag_walk_re << "\t" << stats.diag_walk_im <<"\t NaN \t NaN" << endl;

		myfile.close();
	} catch (ofstream::failure e) {
		cout << "Exception opening/writing file";
	}
}

void QMC_state::Ms_WriteStat_Support(std::string filename, int time)
{
	ofstream myfile;
	stringstream s;

	myfile.exceptions(ofstream::failbit | ofstream::badbit);
	try {

		myfile.open(filename.c_str(), ios::app);

		if(iter%energy_freq == 0)
		{
			myfile << time << "\t" << stats.det_nbr << "\t" << stats.Mx.real() << "\t" << stats.Mx.imag()<< "\t" << stats.My.real() << "\t" << stats.My.imag() << "\t" << stats.Mz.real() << "\t" << stats.Mz.imag() << "\t" << stats.Div.real() << "\t" << stats.Div.imag() << "\t" << stats.diag_det_nbr << "\t" << stats.det_nbr << endl;
		}
		else
			myfile << time << "\t" << stats.det_nbr<< "\t" << "Nan \t NaN \t Nan \t NaN \t Nan \t NaN \t Nan \t Nan" << "\t" << stats.diag_det_nbr << "\t" << stats.det_nbr << endl;

		myfile.close();
	} catch (ofstream::failure e) {
		cout << "Exception opening/writing file";
	}
}

void QMC_state::Ms_Write_Reset(int iter)
{
	ofstream myfile;
	stringstream s;

	s << path << "result/Master_parameter_reset.res";

	myfile.exceptions(ofstream::failbit | ofstream::badbit);
	try {

		myfile.open(s.str().c_str());

		myfile << "Number of processes:" << endl;
		myfile << proc << endl;
		myfile << "Cycle:" << endl;
		myfile << iter << endl;
		myfile << "Is shift calculation turned on?" << endl;
		myfile << pop_ctrl.ctrl << endl;
		myfile << "Actual value of shift:" << endl;
		myfile << pop_ctrl.shift << endl;
		myfile << "Actual number of walkers:" << endl;
		myfile << stats.N_diag << endl;

		myfile.close();
	} catch (ofstream::failure e) {
		cout << "Exception opening/writing file";
	}
}

void QMC_state::Sl_Write_Reset()
{
	ofstream myfile;
	stringstream s;

	s << path << "result/Slave_reset_Proc_" << myrank <<".res";

	myfile.exceptions(ofstream::failbit | ofstream::badbit);
	try {

		myfile.open(s.str().c_str());

		/*for(unordered_map<std::bitset<Bl>, Determinant>::iterator i = determinant.begin(); i != determinant.end(); ++i)
		{
			myfile << IntToBitset(i->second.det) << "\t" << i->second.re << "\t" << i->second.im << endl;
		}*/

		myfile.close();
	} catch (ofstream::failure e) {
		cout << "Exception opening/writing file";
	}
}

/*
void QMC_state::Ms_Read_Reset(string filename, vector<bitset<Bl>> & refdets, vector<int> & walkers) {
	string line;
	stringstream s;
	vector<string> arr;
	ifstream myfile;
	int proc_nbr_read, prc_nbr, walk_nbr_re, walk_nbr_im;
	bitset<Bl> det;

	s << filename << "Master_parameter_reset.res";

	///////////////////////////////////////////////////////////////////////
	//PARAMETERS
	myfile.exceptions(ifstream::failbit | ifstream::badbit);
	try {
		myfile.open(s.str().c_str());

		getline(myfile, line); getline(myfile, line);
		proc_nbr_read = atoi(line.c_str()) - 1;

		getline(myfile, line); getline(myfile, line);
		rest_iter = atoi(line.c_str());

		getline(myfile, line); getline(myfile, line);
		if(atoi(line.c_str()) == 1)
			pop_ctrl.ctrl = 1;

		getline(myfile, line); getline(myfile, line);
		pop_ctrl.shift = atof(line.c_str());

		getline(myfile, line); getline(myfile, line);
		stats.N_diag = atoi(line.c_str());

		myfile.close();

	} catch (ifstream::failure e) {

		cout << "Parameters - Exception opening/reading file" << endl;
	}
	///////////////////////////////////////////////////////////////////////


	///////////////////////////////////////////////////////////////////////
	//DETERMINANTS
	for(int i = 0; i<proc_nbr_read; ++i)
	{
		ifstream myfile;
		s.str("");
		s << filename << "Slave_reset_Proc_" << i << ".res";

		try {
			myfile.open(s.str().c_str());

			while(getline(myfile, line))
			{
				split(line, arr, '\t');

				walk_nbr_re = atoi(arr[1].c_str());
				walk_nbr_im= atoi(arr[2].c_str());
				for(int j = 0; j<Bl; ++j)
					det[Bl - 1 - j] = boost::lexical_cast<int>(arr[0][j]);

				prc_nbr = MPI_Hash(det, proc - 1);
				spawned_send[prc_nbr].push_back({det, walk_nbr_re, walk_nbr_im});
			}

			myfile.close();

		} catch (ifstream::failure e) {

			cout << "Determinants - Exception opening/reading file - Is bitset length correct?" << endl;
		}
	}
	///////////////////////////////////////////////////////////////////////

}*/

/*
void QMC_state::Ms_Restart(MPI_Datatype & pop_c, MPI_Datatype & spawn_type, MPI_Datatype & bits)
{
	vector<bitset<Bl>> refdet;
	vector<int> walk;

	Ms_Read_Reset(restart_path, refdet, walk);

	for(int i=0; i<proc - 1; ++i)
		MPI_Send(&rest_iter, 1, MPI_INT, i, 1, MPI_COMM_WORLD);

	for(int i=0; i<proc - 1; ++i)
		MPI_Send(&stats.N_diag, 1, MPI_INT, i, 2, MPI_COMM_WORLD);

	Ms_SendShift(pop_c);

	for(int i=0; i<proc - 1; ++i)
		MPI_Send(&spawned_send[i][0], spawned_send[i].size(), spawn_type, i, 2, MPI_COMM_WORLD);

}

void QMC_state::Sl_Restart(MPI_Datatype & pop_c, MPI_Datatype & spawn_type, MPI_Datatype & bits, System_Model & sys)
{
	MPI_Status status;

	MPI_Recv(&rest_iter, 1, MPI_INT, master, 1, MPI_COMM_WORLD,	MPI_STATUS_IGNORE);

	MPI_Recv(&stats.N_diag, 1, MPI_INT, master, 2, MPI_COMM_WORLD,	MPI_STATUS_IGNORE);

	MPI_Recv(&pop_ctrl, 1, pop_c, master, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	int length;
	MPI_Probe(master, 2, MPI_COMM_WORLD, &status);
	MPI_Get_count(&status, spawn_type, &length);
	spawned_recv[0].resize(length);
	MPI_Recv(&spawned_recv[0][0], length, spawn_type, master, 2,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	for (unsigned int i = 0; i < spawned_recv[0].size(); ++i) {
		Generate_Determinant(spawned_recv[0][i].det,spawned_recv[0][i].nbr_re, spawned_recv[0][i].nbr_im, sys);
	}
*/

