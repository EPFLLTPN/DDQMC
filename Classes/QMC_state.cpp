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

		myfile.close();

	} catch (ifstream::failure e) {

		cout << "Exception opening/reading file" << endl;
	}

	gsl_rng_env_setup();
	gsl_random = gsl_rng_alloc(gsl_rng_default);
	gsl_rng_set(gsl_random, this->seed);

	//setting up the reference determinant and the starting one -> put it to the correct process
	this->reference = sys.prop.reference;

	ref_process = MPI_Hash(this->reference.to_ulong(), proc-1);

	if(myrank == ref_process)
	{
		Generate_Determinant(this->reference, 100, sys);
	}

	//initialize the temporary spawned walker storage vector
	vector<Spawned_data> temp;
	for(int i=0; i<proc - 1; ++i)
	{
		spawned_send.push_back(temp);
		spawned_recv.push_back(temp);
	}

	//variable initialization
	stats.N_w = 0;
	stats.N_w_old = 1;
	stats.E = 0;
	stats.det_nbr = 0;
	stats.E_div = 0;

	//population control off
	pop_ctrl.ctrl = 0;
	pop_ctrl.alert = 0;
	pop_ctrl.ref_on = 0;

	requests.resize(proc - 1);
}

void QMC_state::Sl_DeterminantStep(System_Model & sys)
{
	Random_excitation rnd_exc;
	int n_sp, proc_nbr, kill;
	double p_death;

	for(int d = 0; d < determinant.size(); ++d)
	{
		for(int w =0; w < abs(walker[d]); ++w)
		{
			//Spawning step
			//Generate random excitation

			rnd_exc = sys.Gen_Rnd_Excitation(gsl_random, determinant[d], det_prop[d].occ_list, det_prop[d].renorm);

			//Try spawning
			n_sp = SpawningAttempt(rnd_exc.H_elm, rnd_exc.p_gen, gsl_random, walker[d]/abs(walker[d]));

			if(n_sp != 0)
			{
					//Store the spawned walkers into the ordered spawn list
				proc_nbr = MPI_Hash(rnd_exc.child, proc - 1);
				spawned_send[proc_nbr].push_back({rnd_exc.child, n_sp});
			}
		}

		//Cloning/killing step
		//dying probability times the whole population of the determinant
		//number that dies for sure plus stochastic part
		p_death = tau * (sys.HamiltonianElement(determinant[d], determinant[d]) - pop_ctrl.shift) * abs((double)walker[d]);

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

		//Find the number of particles
		//If kill is positive, then particles are killed
		if(walker[d] < 0)
			walker[d] += kill;
		else
			walker[d] -= kill;


		stats.N_w += abs(walker[d]);

		if (walker[d] != 0) {

			if (pop_ctrl.ref_on == 0) {
				stats.E += sys.HamiltonianElement(determinant[d], reference)
						* (double) (walker[d]);

				if (determinant[d] == reference)
					stats.E_div += (double) walker[d];

			} else {

				for (int i = 0; i < ref_dets.size(); ++i) {
					stats.E += sys.HamiltonianElement(determinant[d],
							ref_dets[i]) * (double) (walker[d])
							* (double) (ref_walk[i]);

					if (ref_dets[i] == determinant[d])
						stats.E_div += (double) (walker[d])
								* (double) (ref_walk[i]);
				}

			}
		}
	}

	stats.det_nbr = determinant.size();
}

void QMC_state::Sl_SendStat(MPI_Datatype & stat_type)
{
	MPI_Isend(&stats, 1, stat_type, master, 0, MPI_COMM_WORLD, &request);
}

void QMC_state::Ms_RecvStat(MPI_Datatype & stat_type)
{
	Stats stat_temp;

	for(int i = 0; i < proc - 1; ++i)
	{
		MPI_Recv(&stat_temp, 1, stat_type, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		stats.N_w += stat_temp.N_w;
		stats.E += stat_temp.E;
		stats.det_nbr += stat_temp.det_nbr;
		stats.E_div += stat_temp.E_div;
	}

	//energy calculation
	stats.E /= ((double)stats.E_div);

}

void QMC_state::Ms_CalcShift()
{
	//shift calculation
		pop_ctrl.shift = pop_ctrl.shift - (shift_cnst / tau) * log10((double)stats.N_w / (double)stats.N_w_old);

}

void QMC_state::Ms_SendShift(MPI_Datatype & pop_type)
{

	for(int i = 0; i < proc - 1; ++i)
	{
		MPI_Isend(&pop_ctrl, 1, pop_type, i, 4, MPI_COMM_WORLD, &request);
	}

	if(pop_ctrl.alert == 1)
		pop_ctrl.alert = 0;
}

void QMC_state::Sl_RecvShift(MPI_Datatype & pop_type)
{
	MPI_Recv(&pop_ctrl, 1, pop_type, master, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

void QMC_state::Sl_SendSpawned(MPI_Datatype & sp_type)
{
	for(int i=0; i<proc - 1; ++i)
		MPI_Isend(&spawned_send[i][0], spawned_send[i].size(), sp_type, i, 2, MPI_COMM_WORLD, &requests[i]);
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

void QMC_state::Sl_SendRecvRef(MPI_Status & status, System_Model & sys)
{
	//prepare the vector to send the populated determinants
	vector<ulong> send_dets, get_dets;
	vector<int> send_walk, get_walk;

	for(int i = 0; i<determinant.size(); ++i)
	{
		send_dets.push_back(determinant[i].to_ulong());
		send_walk.push_back(walker[i]);
	}

	for(int i=0; i<proc - 1; ++i)
	{
		MPI_Isend(&send_dets[0], send_dets.size(), MPI_UNSIGNED_LONG, i, 7, MPI_COMM_WORLD, &ref1);
		MPI_Isend(&send_walk[0], send_walk.size(), MPI_INT, i, 8, MPI_COMM_WORLD, &ref2);
	}

	int length;
	for(int i = 0; i < proc - 1; ++i)
	{
		MPI_Probe(MPI_ANY_SOURCE, 7, MPI_COMM_WORLD, &status);
		MPI_Get_count(&status, MPI_UNSIGNED_LONG, &length);
		get_dets.resize(length);
		MPI_Recv(&get_dets[0], length, MPI_UNSIGNED_LONG, status.MPI_SOURCE, 7, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		for(int j=0; j<get_dets.size();++j)
		{
			bitset<B> de(get_dets[j]);
			ref_dets.push_back(de);
		}

		get_walk.resize(length);
		MPI_Recv(&get_walk[0], length, MPI_INT,	status.MPI_SOURCE, 8, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		for(int j=0; j<get_walk.size();++j)
		{
			ref_walk.push_back(get_walk[j]);
		}
	}
}

void QMC_state::Sl_SortSpawned(int proc, System_Model & sys)
{
	for(int i = 0; i < spawned_recv[proc].size(); ++i)
	{
		bitset<B> dd(spawned_recv[proc][i].det);
		find_det = location.find(dd);

		if(find_det == location.end())
		{
			Generate_Determinant(dd, spawned_recv[proc][i].nbr, sys);
		}
		else
		{
			walker[find_det->second] += spawned_recv[proc][i].nbr;
		}
	}
}




void QMC_state::Generate_Determinant(bitset<B> det, int nbr_w, System_Model & sys)
{
	determinant.push_back(det);
	walker.push_back(nbr_w);
	location.insert({det, determinant.size()-1});

	Determinant_properties properties;
	properties.renorm = sys.Calc_Occ_List_Renorm(det, properties.occ_list);

	det_prop.push_back(properties);
}

int QMC_state::SpawningAttempt(double hamelm, double pgen, gsl_rng * rand_gen, int parent_sign)
{
	double p_spawn;
	int n_spawn;

	//Calculate the probability that spawning is successful
	p_spawn = tau * abs(hamelm) / pgen;

	//Need to take into account that multiple spawn can appear
	n_spawn = floor(p_spawn);
	p_spawn -= n_spawn;


	if(p_spawn > gsl_rng_uniform(rand_gen))
		++n_spawn;

	//decide the sign of the spawned walkers
	if(n_spawn > 0)
	{
		//If H_ij is positive, then the spawned walker is of opposite sign to the parent
		// otherwise the spawned walked is of the same sign as the parent
		if(hamelm > 0)
			n_spawn = n_spawn * parent_sign * (-1);
		else
			n_spawn = n_spawn * parent_sign;
	}

	return n_spawn;
}

void QMC_state::Cleaner()
{
	stats.N_w_old = stats.N_w;
	stats.N_w = 0;
	stats.E = 0;
	stats.det_nbr = 0;
	stats.E_div = 0;

	for (int i = 0; i < proc - 1; ++i)
	{
		spawned_send[i].clear();
		spawned_recv[i].clear();
	}
}

void QMC_state::Ms_WriteStat(string filename, int time)
{
	ofstream myfile;
	stringstream s;

	myfile.exceptions(ofstream::failbit | ofstream::badbit);
	try {

		myfile.open(filename.c_str(), ios::app);

		myfile << time <<" \t " << stats.N_w <<" \t " << pop_ctrl.shift <<" \t " << stats.E << " \t " << stats.det_nbr << endl;


		myfile.close();
	} catch (ofstream::failure e) {
		cout << "Exception opening/writing file";
	}
}

