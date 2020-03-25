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
#include  "Header/Lindbladian_Test.h"

using namespace std;

//Global variables
int proc, myrank, master, energy_freq;
MPI_Request request;
int rest_iter = 0;
string path;
string restart_path = "";
string res_path;
bool restart = false;
int iter;
complex<double> im;
int InitLimit;

//Functions
string ModelType(string filename);
void FinalTime(string filename, clock_t & begin);
void DataTest(int nbr, int walker, int det);

int main(int argc, char *argv[]) {

	complex<double> imhelp(0,1);
		im = imhelp;
		//initialize variable
		System_Model sys;
		QMC_state * qmc_state;
		MPI_Datatype properties, stat_type, spawn_type, pop_ctrl, bits_type;
		//Statistics statistic;
		string outPath, outPut, outPut_Sup;
		clock_t begin_time = clock();
		int mainrun = 1;

		//support variables
		stringstream s;
		string syssize, yesno, dir;
		string JJ, ext, iii;

		path = argv[1];
		res_path = argv[2];
		syssize = argv[3];
		JJ = argv[4];
		yesno = argv[5];
		dir = argv[6];
		ext = argv[7];
		iii = argv[8];


		s.str("");
		if( yesno == "NO")
			s << path << "/model/model_" << syssize << "_" << JJ << "_NO_" << iii << ".txt";
		else
			s << path << "/model/model_" << syssize << "_" << JJ << "_" << dir << "_" << ext << "_" << iii << ".txt";

	//if(restart_path != "-")
	//	restart = true;

	//initialize MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &proc);
	MPI_Status status;
	master = proc - 1;

	//initialize Data type
	sys.MPI_Data(properties, stat_type, spawn_type, pop_ctrl, bits_type);

	if (myrank == master)
		sys.Init(s.str());

	s.str("");
		if( yesno == "NO")
			s << path << "/config/config_" << syssize << "_" << JJ << "_NO_" << iii << ".txt";
		else
			s << path << "/config/config_" << syssize << "_" << JJ << "_" << dir << "_" << ext << "_" << iii << ".txt";


	MPI_Bcast(&sys.prop, 1, properties, master, MPI_COMM_WORLD);
	//Initializing the qmc_state object (everyone, even master will need the config data)

	qmc_state = new QMC_state(s.str(), sys);

	//Write the first data into the output files by MASTER
	if(myrank == master){
		outPut = sys.OutputGen(qmc_state->cycle, qmc_state->target_nbr, qmc_state->pop_ctrl.shift, qmc_state->reference, qmc_state->tau, outPath, qmc_state->seed, qmc_state->shift_cnst, qmc_state->impS);
		outPut_Sup = outPut;
		outPut_Sup.erase(outPut_Sup.size()-4);
		s.str("");
		s << outPut_Sup << "_SUPPORT.out";
		outPut_Sup = s.str();


		if(mainrun == 0)
		{
			/*
			//Lindbladian Test generator
			Lindbladian_Test * test = new Lindbladian_Test(sys, qmc_state->tau);

			for (int timestep = 0; timestep < 10000; ++timestep){
				test->Generate_Density_Matrix();
			}
			test->PrintDM();
			*/

		/*	bitset<Bl> proba(string("011001"));
			qmc_state->Generate_Determinant(proba, 1000, 1000, sys);

			qmc_state->TestPerformance(qmc_state->determinant.find(proba)->second, 1);
			for(unsigned int walk=500; walk<=50000;walk+=500)
				qmc_state->TestPerformance(qmc_state->determinant.find(proba)->second, walk);
				*/
		}
	}

	if(mainrun == 1)
	{
	////////////////////////////////////////////////////////////////////////////////
	//RESTART
	if(restart)
	{
		if(myrank == master)
		{
			//qmc_state->Ms_Restart(pop_ctrl, spawn_type, bits_type);
		}
		else
		{
			//qmc_state->Sl_Restart(pop_ctrl, spawn_type, bits_type, sys);
		}

		qmc_state->Cleaner();
	}


	////////////////////////////////////////////////////////////////////////////////
	//MAIN CYCLE

	for(iter = rest_iter; iter < qmc_state->cycle; ++iter)	{

		//We have to distinguish between master and slave behaviour
		if(myrank == master)
		{
			qmc_state->Ms_RecvStat(stat_type);

			if(qmc_state->pop_ctrl.ctrl == 0)
			{
				if( (qmc_state->stats.N_diag) >= qmc_state->target_nbr)
				{
					qmc_state->pop_ctrl.ctrl = 1;
				}
			}
			else
			{
				qmc_state->Ms_CalcShift();
			}

			qmc_state->Ms_SendShift(pop_ctrl);
			qmc_state->Ms_WriteStat(outPut, iter);
			qmc_state->Ms_WriteStat_Support(outPut_Sup, iter);

			//if(iter % qmc_state->write_out_period == 0)
			//	qmc_state->Ms_Write_Reset(iter);

			if(iter%energy_freq ==0)
			//	cout << iter <<" \t " << qmc_state->stats.N_tot <<" \t " <<  qmc_state->stats.det_nbr << "\t" <<
				cout << iter <<" \t " << qmc_state->stats.N_tot <<" \t " << qmc_state->stats.N_diag << "\t"
				<<  qmc_state->pop_ctrl.shift << "\t" << qmc_state->stats.det_nbr <<" \t "
				 << qmc_state->stats.Mx/ qmc_state->stats.Div << "\t "
				<< qmc_state->stats.My/ qmc_state->stats.Div   << "\t"
				<< qmc_state->stats.Mz/ qmc_state->stats.Div   << "\t" <<
			//	getCurrentRSS() << "      " << getPeakRSS() << endl;
				qmc_state->stats.diag_walk_re << "\t" << qmc_state->stats.diag_walk_im << endl;
			else
				cout << iter <<" \t " << qmc_state->stats.N_tot << "\t"<< qmc_state->stats.N_diag << "\t " <<  qmc_state->pop_ctrl.shift << "\t"<< qmc_state->stats.det_nbr << "\t - \t - \t - \t"<< qmc_state->stats.diag_walk_re << "\t" << qmc_state->stats.diag_walk_im << endl;

		}
		else
		{
			qmc_state->Sl_DeterminantStep(sys);
			qmc_state->Sl_SendStat(stat_type); //non-blocking
			qmc_state->Sl_SendSpawned(spawn_type);  //non-blocking
			qmc_state->Sl_RecvSpawned(spawn_type,status,sys);  //blocking
			qmc_state->Sl_RecvShift(pop_ctrl); //blocking

			/*if(iter % qmc_state->write_out_period == 0)
			{
				qmc_state->Sl_Write_Reset();
			}*/


			//MPI_Waitall(proc - 1, qmc_state->ReqArray, MPI_STATUSES_IGNORE);

			//DataTest(myrank, qmc_state->stats.N_diag, qmc_state->determinant.size());

		}
		MPI_Barrier(MPI_COMM_WORLD);


		qmc_state->Cleaner();
	}

	if(myrank == master)
	{
		complex<double> temp(0,0), trace(0,0);

		for(int i=0; i<proc - 1; ++i)
		{
			MPI_Recv(&temp, 1, MPI_CXX_DOUBLE_COMPLEX, i, 11, MPI_COMM_WORLD,	MPI_STATUS_IGNORE);
			trace += temp;
		}

		for(int i=0; i<proc - 1; ++i){
			MPI_Send(&trace, 1, MPI_CXX_DOUBLE_COMPLEX, i, 12, MPI_COMM_WORLD);
		}
	}
	else
	{
		complex<double> tr_count(0,0), trace;

		unordered_map<RepType, Determinant, ArrayHasher>::iterator d = qmc_state->determinant.begin();
		while(d != qmc_state->determinant.end())
		{
			if(d->second.diagonal)
			{
				complex<double> elm(d->second.re,  d->second.im);
				tr_count += elm;
			}

			++d;
		}

		MPI_Send(&tr_count, 1, MPI_CXX_DOUBLE_COMPLEX, master, 11, MPI_COMM_WORLD);
		MPI_Recv(&trace, 1, MPI_CXX_DOUBLE_COMPLEX, master, 12, MPI_COMM_WORLD,	MPI_STATUS_IGNORE);

		qmc_state->PrintDensity(trace);
	}


	if(myrank == master)
		FinalTime(outPut, begin_time);

	}

	sys.MPI_Free(properties,stat_type,spawn_type, pop_ctrl, bits_type);
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

void FinalTime(string filename, clock_t & begin)
{
	ofstream myfile;
	stringstream s;

	time_t now = time(0);

	myfile.exceptions(ofstream::failbit | ofstream::badbit);
	try {

		myfile.open(filename.c_str(), ios::app);

		myfile << "Finishing time:\t" << ctime(&now) <<endl;
		myfile << "Execution time: \t" << (float( clock () - begin ) /  CLOCKS_PER_SEC)/60.0 << "  hours" << endl;

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

