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
 * QMC_state.cpp
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

		myfile.close();

	} catch (ifstream::failure e) {

		cout << "Exception opening/reading file" << endl;
	}

	gsl_rng_env_setup();
	gsl_random = gsl_rng_alloc(gsl_rng_default);
	gsl_rng_set(gsl_random, this->seed);

	//setting up the reference determinant and the starting one -> put it to the correct process
	this->reference = sys.prop.reference;

	ref_process = MPI_Hash(this->reference, proc-1);

	if(myrank == ref_process)
	{
		if(!restart)
			Generate_Determinant(this->reference, 150, 0, sys);
	}

	//initialize the temporary spawned walker storage vector
	spawned_send.resize(proc - 1);
	spawned_recv.resize(proc - 1);

	//variable initialization
	stats.N_w = 0;
	stats.N_w_old = 150;
	stats.N_tot = 0;
	stats.det_nbr = 0;
	stats.diag_walk_re = 0;
	stats.diag_walk_im = 0;
	//stats.d_r = 0;
	//stats.d_i = 0;
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
	//int dcr = 0, dci = 0;
	double parent_sign_re, parent_sign_im;
	int n_sp;
	unsigned int * sp;


	unordered_map<std::bitset<Bl>, Determinant>::iterator d = determinant.begin();
	while(d != determinant.end())
	{

		// REAL PARENT
		if(abs(d->second.re) != 0)
		{
			parent_sign_re = d->second.re / abs(d->second.re);

			if (d->second.prop.P_tot > 1)
				n_sp = floor(d->second.prop.P_tot)+ gsl_ran_binomial(gsl_random,d->second.prop.P_tot - floor(d->second.prop.P_tot), abs(d->second.re));
			else
				n_sp = gsl_ran_binomial(gsl_random, d->second.prop.P_tot, abs(d->second.re));

			//MULTINOMIAL
			if (n_sp > 0)
			{
				sp = new unsigned int[d->second.prop.conn_prob.size()];
				gsl_ran_multinomial(gsl_random, d->second.prop.conn_prob.size(), n_sp,	&d->second.prop.conn_prob[0], sp);

				for (unsigned int i = 0; i < d->second.prop.n_im; ++i)
					if (sp[i] != 0)
					{
						proc_nbr = MPI_Hash(d->second.prop.conn_states[i], proc - 1);
						spawned_send[proc_nbr].push_back({ d->second.prop.conn_states[i], 0,parent_sign_re*((d->second.prop.conn_ampl[i].imag())/(abs(d->second.prop.conn_ampl[i].imag())))* sp[i] });
					}

				for (unsigned int i = d->second.prop.n_im;i < d->second.prop.conn_states.size(); ++i)
					if (sp[i] != 0)
					{
						proc_nbr = MPI_Hash(d->second.prop.conn_states[i], proc - 1);
						spawned_send[proc_nbr].push_back({ d->second.prop.conn_states[i],parent_sign_re*((d->second.prop.conn_ampl[i].real())/ (abs(d->second.prop.conn_ampl[i].real())))* sp[i], 0 });
					}
			}
		}

		// IMAGINARY PARENT
		if(abs(d->second.im) != 0)
		{
			parent_sign_im = d->second.im / abs(d->second.im);

			if (d->second.prop.P_tot > 1)
				n_sp = floor(d->second.prop.P_tot)+ gsl_ran_binomial(gsl_random,d->second.prop.P_tot - floor(d->second.prop.P_tot), abs(d->second.im));
			else
				n_sp = gsl_ran_binomial(gsl_random, d->second.prop.P_tot, abs(d->second.im));

			//MULTINOMIAL
			if (n_sp > 0)
			{
				sp = new unsigned int[d->second.prop.conn_prob.size()];
				gsl_ran_multinomial(gsl_random, d->second.prop.conn_prob.size(), n_sp,	&d->second.prop.conn_prob[0], sp);

				for (unsigned int i = 0; i < d->second.prop.n_im; ++i)
					if (sp[i] != 0)
					{
						proc_nbr = MPI_Hash(d->second.prop.conn_states[i], proc - 1);
						spawned_send[proc_nbr].push_back({ d->second.prop.conn_states[i], parent_sign_im*(-1.0)*((d->second.prop.conn_ampl[i].imag())/(abs(d->second.prop.conn_ampl[i].imag())))* sp[i],0 });
					}

				for (unsigned int i = d->second.prop.n_im;i < d->second.prop.conn_states.size(); ++i)
					if (sp[i] != 0)
					{
						proc_nbr = MPI_Hash(d->second.prop.conn_states[i], proc - 1);
						spawned_send[proc_nbr].push_back({ d->second.prop.conn_states[i],0,parent_sign_im*((d->second.prop.conn_ampl[i].real())/ (abs(d->second.prop.conn_ampl[i].real())))* sp[i] });
					}
			}
		}

		//Cloning/killing step
		//dying probability times the whole population of the determinant
		//number that dies for sure plus stochastic part
		complex<double> d_a;
		if(d->second.diagonal){
			d_a = d->second.prop.diag_ampl - pop_ctrl.shift;
			stats.N_w += abs(d->second.re) + abs(d->second.im);
		}
		else
			d_a = d->second.prop.diag_ampl;
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
				//dcr += d->second.re;
				//dci += d->second.im;
				diag_count_re += abs(d->second.re);
				diag_count_im += abs(d->second.im);}

			++d;
		}
		else
		{
			determinant.erase(d++);
		}
	}

	//stats.d_r = dcr;
	//stats.d_i = dci;
	stats.diag_walk_re = diag_count_re;
	stats.diag_walk_im = diag_count_im;
	stats.det_nbr = determinant.size();
}

void QMC_state::Sl_SendStat(MPI_Datatype & stat_type)
{

	MPI_Isend(&stats, 1, stat_type, master, 0, MPI_COMM_WORLD, &request);
}

void QMC_state::PrintDensity(complex<double> trace)
{
	unordered_map<std::bitset<Bl>, Determinant>::iterator d = determinant.begin();
	while(d != determinant.end())
	{
		if(d->second.diagonal == 1)
		{
			complex<double> elm(d->second.re, d->second.im);
			cout << d->second.det <<"   "<< (elm/abs(trace)).real()  << "    " << (elm/abs(trace)).imag() << endl;
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
		stats.N_w += stat_temp.N_w;
		stats.N_tot += stat_temp.N_tot;
		stats.det_nbr += stat_temp.det_nbr;
		stats.diag_walk_re += stat_temp.diag_walk_re;
		stats.diag_walk_im += stat_temp.diag_walk_im;
		stats.Mx += stat_temp.Mx;
		stats.My += stat_temp.My;
	//	stats.d_r += stat_temp.d_r;
		//stats.d_i += stat_temp.d_i;
		stats.Mz += stat_temp.Mz;
		stats.Div += stat_temp.Div;
	}

}
int QMC_state::SpawningAttempt(complex<double> & ampl, double pgen, gsl_rng * rand_gen, int parent, int child, int par_sign)
{
	double p_spawn;
	int n_spawn;
	int sign;

	//Calculate the probability that spawning is successful
	if(parent) //real parent
	{
		if(child){ //trying to spawn real child
			p_spawn = tau * abs(ampl.real()) / pgen;
			sign =  (ampl.real())/abs(ampl.real()) * par_sign;}
		else{   //trying to spawn imaginary child
			p_spawn = tau * abs(ampl.imag()) / pgen;
			sign =  (ampl.imag())/abs(ampl.imag()) * par_sign;}
	}
	else   //imaginary parent
	{
		if(child){ //trying to spawn real child
			p_spawn = tau * abs(ampl.imag()) / pgen;
			sign =  (-1.0)*(ampl.imag())/abs(ampl.imag()) * par_sign;}
		else{   //trying to spawn imaginary child
			p_spawn = tau * abs(ampl.real()) / pgen;
			sign =  (ampl.real())/abs(ampl.real()) * par_sign;}
	}

	//Need to take into account that multiple spawn can appear
	n_spawn = floor(p_spawn);
	p_spawn -= n_spawn;

	if(p_spawn > gsl_rng_uniform(rand_gen))
		++n_spawn;

	//decide the sign of the spawned walkers
	if(n_spawn > 0)
		n_spawn = n_spawn * sign;

	return n_spawn;
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
		kill = floor(p_death)+ gsl_ran_binomial(gsl_random,p_death - floor(p_death),walker);
	} else
		kill = gsl_ran_binomial(gsl_random, p_death,walker);

	if (kill > 0)
		kill = kill * sign;

	return kill;
}

void QMC_state::Ms_CalcShift()
{
	//shift calculation
		pop_ctrl.shift = pop_ctrl.shift + (shift_cnst / tau) * log10((double)stats.N_w / (double)stats.N_w_old);

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
			Generate_Determinant(spawned_recv[proc][i].det, spawned_recv[proc][i].nbr_re, spawned_recv[proc][i].nbr_im, sys);
		}
		else
		{
			find_det->second.re += spawned_recv[proc][i].nbr_re;
			find_det->second.im += spawned_recv[proc][i].nbr_im;
		}
	}
}


void QMC_state::Generate_Determinant(bitset<Bl> det, int nbr_re, int nbr_im, System_Model & sys)
{
	int count;
	vector<int> where;
	Determinant d;

	d.det = det;
	d.re = nbr_re;
	d.im = nbr_im;

	Determinant_properties properties;
	properties.P_tot = 0;
	properties.n_im = 0;

	////////////////////////////////////////////////////////////////////////////////
	//Calculating estimator element for the given density matrix state
	count = 0;
	for(int N = 0; N<sys.prop.xSite;++N)
		if(det[N]^det[N+sys.prop.xSite]){
			++count;
			where.push_back(N);}

	//DIAGONAL
	if(count == 0)
	{
		d.diagonal = 1;
		for(int i = 0; i<sys.prop.xSite; ++i)
			properties.mz += 0.5 * pow(-1.0, det[i] + 1);

		properties.mz /= (double)sys.prop.xSite;
	}
	//OFF-DIAGONAL
	else if(count == 1)
	{
		d.diagonal = 0;

		properties.mx = 0.5;

		if(det[where[0]] == 0)
			properties.my = -0.5*im;
		else
			properties.my = 0.5*im;

		properties.my /= (double)sys.prop.xSite;
	}
	else
		d.diagonal = 0;

	////////////////////////////////////////////////////////////////////////////////
	//Filling up all the possible connections
	bitset<Bl> state_l;
	int left;
	int help_init;

	//COL
	for(int init = 0; init<sys.prop.xSite;++init)
	{
		state_l = det;

		help_init = init;
		left = (sys.prop.xSite + (help_init - 1) % sys.prop.xSite) % sys.prop.xSite;

		//////////////////////////////////////////////////////////////////////
		//Hopping type + Double hopping type
		state_l.flip(left);
		state_l.flip(init);

		if(det[left] != det[init]){
			properties.conn_states.push_back(state_l);
			properties.n_im += 1;

			properties.conn_ampl.push_back(im*sys.prop.hop);
			properties.conn_prob.push_back(abs(im*sys.prop.hop)*tau);
			properties.P_tot += abs(im*sys.prop.hop)*tau;
		}
		else{
			if(sys.prop.double_hop != 0)
			{
				properties.conn_states.push_back(state_l);
				properties.n_im += 1;

				properties.conn_ampl.push_back(im*sys.prop.double_hop);
				properties.conn_prob.push_back(abs(im*sys.prop.double_hop)*tau);
				properties.P_tot += abs(im*sys.prop.double_hop)*tau;
			}
		}
	}

	//ROW
	for(int init = sys.prop.xSite; init<2*sys.prop.xSite;++init)
	{
		state_l = det;

		help_init = init-sys.prop.xSite;
		left = (sys.prop.xSite + (help_init - 1) % sys.prop.xSite) % sys.prop.xSite;
		left += sys.prop.xSite;

		//////////////////////////////////////////////////////////////////////
		//Hopping type + Double hopping type
		state_l.flip(left);
		state_l.flip(init);

		if(det[left] != det[init]){
			properties.conn_states.push_back(state_l);
			properties.n_im += 1;

			properties.conn_ampl.push_back(-im*sys.prop.hop);
			properties.conn_prob.push_back(abs(-im*sys.prop.hop)*tau);
			properties.P_tot += abs(-im*sys.prop.hop)*tau;
		}
		else{
			if(sys.prop.double_hop != 0)
			{
				properties.conn_states.push_back(state_l);
				properties.n_im += 1;

				properties.conn_ampl.push_back(-im*sys.prop.double_hop);
				properties.conn_prob.push_back(abs(-im*sys.prop.double_hop)*tau);
				properties.P_tot += abs(-im*sys.prop.double_hop)*tau;
			}
		}
	}

	//ANY
    for(int init = 0; init < sys.prop.xSite; ++init)
    {
    	if((det[init] == 1) && (det[init + sys.prop.xSite] == 1))
    	{
    	  	state_l = det;
    	  	state_l.flip(init);
    	  	state_l.flip(init + sys.prop.xSite);

    	  	properties.conn_states.push_back(state_l);
    	  	properties.conn_ampl.push_back(sys.prop.gamma);
    	  	properties.conn_prob.push_back(abs(sys.prop.gamma)*tau);
    	  	properties.P_tot += abs(sys.prop.gamma)*tau;
    	}
    }

    for(unsigned int i = 0; i<properties.conn_prob.size(); ++i)
    	properties.conn_prob[i] /= properties.P_tot;

    double h_ii = sys.DiagonalElement(det, 0);
    double h_jj = sys.DiagonalElement(det, 1);
    properties.diag_ampl = im*(h_jj-h_ii)-0.5 * sys.prop.gamma * det.count();

	d.prop = properties;

/*
	cout << d.det;
	cout << "       All" << endl;
	for(int i = 0; i<d.prop.conn_states.size();++i)
		cout <<"        "<<d.prop.conn_states[i] <<"        " << d.prop.conn_ampl[i] <<"         "<< d.prop.diag_ampl<< endl;
*/


	determinant.insert({det, d});
}



void QMC_state::Cleaner()
{
	stats.N_w_old = stats.N_w;
	stats.N_w = 0;
	stats.N_tot = 0;
	stats.det_nbr = 0;
	//stats.d_r = 0;
	//stats.d_i = 0;
	stats.diag_walk_re = 0;
	stats.diag_walk_im = 0;
	stats.Div = 0;
	stats.Mx = 0;
	stats.My = 0;
	stats.Mz = 0;

	for (int i = 0; i < proc - 1; ++i)
	{
		spawned_send[i].clear();
		spawned_recv[i].clear();
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
			myfile << time <<" \t " << stats.N_tot <<" \t " << stats.N_w << "\t" <<
			pop_ctrl.shift <<" \t " << stats.det_nbr << "\t" <<
			stats.Mx / stats.Div  << "\t" <<
			stats.My / stats.Div  << "\t" <<
			stats.Mz / stats.Div << "\t" <<
			stats.diag_walk_re << "\t" << stats.diag_walk_im << "\t" <<
			(stats.Mz/stats.Div).real() << "\t" << (stats.Mz/stats.Div).imag() << endl;
		}
		else
			myfile << time <<" \t " << stats.N_w << "\t"<< stats.diag_walk_re + stats.diag_walk_im << "\t" << pop_ctrl.shift <<" \t " << stats.det_nbr << "\t NaN \t NaN \t NaN \t"<< stats.diag_walk_re << "\t" << stats.diag_walk_im <<"\t NaN \t NaN" << endl;

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
			myfile << time <<" \t " << stats.N_tot <<" \t " << stats.N_w << "\t"<<
			pop_ctrl.shift <<" \t "
			<< stats.det_nbr << "\t" <<
			stats.Mx / stats.Div  << "\t" <<
			stats.My / stats.Div  << "\t" <<
			stats.Mz / stats.Div << "\t" <<
			stats.diag_walk_re << "\t" << stats.diag_walk_im << "\t" <<
			(stats.Mz/stats.Div).real() << "\t" << (stats.Mz/stats.Div).imag() << endl;
		}
		else
			myfile << time <<" \t " << stats.N_w << "\t"<< stats.diag_walk_re + stats.diag_walk_im << "\t" << pop_ctrl.shift <<" \t " << stats.det_nbr << "\t NaN \t NaN \t NaN \t NaN \t NaN \t NaN \t NaN \t"<<stats.diag_walk_re << "\t" << stats.diag_walk_im<<"\t NaN \t NaN" << endl;

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
		myfile << stats.N_w << endl;

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

		for(unordered_map<std::bitset<Bl>, Determinant>::iterator i = determinant.begin(); i != determinant.end(); ++i)
		{
			myfile << i->second.det << "\t" << i->second.re << "\t" << i->second.im << endl;
		}

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
		stats.N_w = atoi(line.c_str());

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
		MPI_Send(&stats.N_w, 1, MPI_INT, i, 2, MPI_COMM_WORLD);

	Ms_SendShift(pop_c);

	for(int i=0; i<proc - 1; ++i)
		MPI_Send(&spawned_send[i][0], spawned_send[i].size(), spawn_type, i, 2, MPI_COMM_WORLD);

}

void QMC_state::Sl_Restart(MPI_Datatype & pop_c, MPI_Datatype & spawn_type, MPI_Datatype & bits, System_Model & sys)
{
	MPI_Status status;

	MPI_Recv(&rest_iter, 1, MPI_INT, master, 1, MPI_COMM_WORLD,	MPI_STATUS_IGNORE);

	MPI_Recv(&stats.N_w, 1, MPI_INT, master, 2, MPI_COMM_WORLD,	MPI_STATUS_IGNORE);

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

