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
 * QMC_main.cpp
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
#include "Header/QMC_state.h"



//Global variables
int proc, myrank, master;
MPI_Request request;
string path;

//Functions
string ModelType(string filename);
void FinalTime(string filename);
void DataTest(int nbr, int walker, int det);

int main(int argc, char *argv[]) {

	//initialize variable
	System_Model sys;
	QMC_state * qmc_state;
	MPI_Datatype properties, stat_type, spawn_type, pop_ctrl;
	Statistics statistic;
	string outPath, outPut;

	//bitstring size
	#undef B
	#define B atoi(argv[2])

	//support variables
	stringstream s;
	path = argv[1];
	s << path << "model.txt";
	double analytical;

	//initialize MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &proc);
	MPI_Status status;
	master = proc - 1;

	//initialize Data type
	sys.MPI_Data(properties, stat_type, spawn_type, pop_ctrl);

	if (myrank == master)
	{
		sys.Init(s.str());
		analytical = sys.Analytical(16);
	}
	MPI_Bcast(&sys.prop, 1, properties, master, MPI_COMM_WORLD);

	//Initializing the qmc_state object (everyone, even master will need the config data)
	s.str(string()); s << path << "config.txt";
	qmc_state = new QMC_state(s.str(), sys);


	//Write the first data into the output files by MASTER
	if(myrank == master){
		outPut = sys.OutputGen(qmc_state->cycle, qmc_state->target_nbr, qmc_state->pop_ctrl.shift, qmc_state->reference, qmc_state->tau, outPath, qmc_state->seed, qmc_state->shift_cnst);
		//statistic.Write_Balance_init(outPath);
	}

	////////////////////////////////////////////////////////////////////////////////
	//MAIN CYCLE
	for(int iter = 0; iter < qmc_state->cycle; ++iter)	{
		//We have to distinguish between master and slave behaviour
		if(myrank == master)
		{
			qmc_state->Ms_RecvStat(stat_type);

			if(qmc_state->pop_ctrl.ctrl == 0)
			{
				if(qmc_state->stats.N_w >= qmc_state->target_nbr)
				{
					qmc_state->pop_ctrl.ctrl = 1;
				}
			}
			else
			{
				qmc_state->Ms_CalcShift();
			}

			//turning on ref_det
			if(qmc_state->stats.det_nbr >= 100)
			{
				qmc_state->pop_ctrl.alert = 1;
				qmc_state->pop_ctrl.ref_on = 1;
			}

			qmc_state->Ms_SendShift(pop_ctrl);
			qmc_state->Ms_WriteStat(outPut, iter);

			cout << iter <<" \t " << qmc_state->stats.N_w <<" \t " << qmc_state->pop_ctrl.shift <<" \t " << qmc_state->stats.E << "\t" << qmc_state->stats.det_nbr << endl;

		}
		else
		{
			qmc_state->Sl_DeterminantStep(sys);
			qmc_state->Sl_SendStat(stat_type); //non-blocking
			qmc_state->Sl_SendSpawned(spawn_type);  //non-blocking
			qmc_state->Sl_RecvSpawned(spawn_type,status,sys);  //blocking
			qmc_state->Sl_RecvShift(pop_ctrl); //blocking

			if(qmc_state->pop_ctrl.alert == 1)
				qmc_state->Sl_SendRecvRef(status, sys);


		}

		if(myrank != master)
		{
			MPI_Request arr [proc-1];
			copy(qmc_state->requests.begin(), qmc_state->requests.end(), arr);
			MPI_Waitall(proc-1, arr, MPI_STATUSES_IGNORE);
			if (qmc_state->pop_ctrl.alert == 1) {
				MPI_Wait(&qmc_state->ref1, MPI_STATUS_IGNORE);
				MPI_Wait(&qmc_state->ref2, MPI_STATUS_IGNORE);
			}

			//DataTest(myrank, qmc_state->stats.N_w, qmc_state->determinant.size());
		}

		qmc_state->Cleaner();
	}


	if(myrank == master)
		FinalTime(outPut);

	sys.MPI_Free(properties,stat_type,spawn_type, pop_ctrl);
	MPI_Finalize();

	return 0;
}

void DataTest(int nbr, int walker, int det)
{
	ofstream myfile;
	stringstream s;

	s << path << "result/" << "Proc_" << nbr <<".out";

	myfile.exceptions(ofstream::failbit | ofstream::badbit);
	try {

		myfile.open(s.str().c_str(), ios::app);

		myfile << walker << "\t" << det << endl;

		myfile.close();
	} catch (ofstream::failure e) {
		cout << "Exception opening/writing file";
	}
}

void FinalTime(string filename)
{
	ofstream myfile;
	stringstream s;

	time_t now = time(0);

	myfile.exceptions(ofstream::failbit | ofstream::badbit);
	try {

		myfile.open(filename.c_str(), ios::app);

		myfile << "Finishing time:\t" << ctime(&now) <<endl;

		myfile.close();
	} catch (ofstream::failure e) {
		cout << "Exception opening/writing file";
	}
}


string ModelType(string filename) {
	string line;
	ifstream myfile;
	myfile.exceptions(ifstream::failbit | ifstream::badbit);
	try {
		myfile.open(filename.c_str());

		getline(myfile, line);
		getline(myfile, line);
		getline(myfile, line);

		myfile.close();
	} catch (ifstream::failure e) {
		cout << "Exception opening/reading file";
	}

	return line;
}

