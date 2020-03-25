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
 * QMC_state.h
 *
 *  Created on: Jul 5, 2016
 *      Author: Alexandra Nagy
 *              alex.nagy.phys@gmail.com
 *
 *      REFERENCES:
 *      - A. Nagy, V. Savona, Driven-dissipative quantum Monte Carlo method for open quantum systems, Physical Review A 97 (5), 052129, 2018
 *      - A. Nagy's PhD Dissertation, École polytechnique fédérale de Lausanne, Quantum Monte Carlo Approach to the non-equilibrium steady state of open quantum systems, 2020
 */

//For every different type of model, one has to rewrite System_Model.h and System_Model.cpp
#include "System_Model.h"
//#include "Main.h"

class QMC_state
{
public:

	QMC_state();
	QMC_state(std::string filename, System_Model & sys);

	gsl_rng * gsl_random;
	std::unordered_map<RepType, Determinant, ArrayHasher>::iterator find_det;
	MPI_Request * ReqArray;

	double tau; //times step
	double impS; // imp_sampl_coeff
	Pop_control pop_ctrl; //boolean for controlling the population and energy shift to be updated with
	double target_nbr; //target particle number
	double shift_cnst; //the constant used for the shift calculation
	std::bitset<Bl> reference; //the reference determinant
	int ref_process; //the process which will contain the reference determinant
	Stats stats; //nbr of walkers on the process and energy estimator to be send to the master
	int cycle; //number of iterations
	int seed; //seed for the random generator
	int write_out_period; //how often it gotta write the result out

	std::unordered_map<RepType, Determinant, ArrayHasher> determinant; //list of the flags of being an initiator

	std::vector<std::vector<Spawned_data>> spawned_send; //vector of vectors associated to the correct process to be sent
	std::vector<std::vector<Spawned_data>> spawned_recv; //vector of vectors associated to the correct process to be received

	void Sl_DeterminantStep(System_Model & sys); //loop through determinants, kill, clone, spawn into temporary

	void Sl_SendStat(MPI_Datatype & stat_type); //send stat. to master (N_walk, E, timer from earlier) - non-blocking
	void Sl_SendSpawned(MPI_Datatype & sp_type); //sort then send spawned to slaves - non-blocking
	void Sl_RecvSpawned(MPI_Datatype & sp_type, MPI_Status & status, System_Model & sys); //receive and organize spawned - blocking
	void Sl_SortSpawned(int proc, System_Model & sys);   //sort out the received spawned walkers
	void Sl_RecvShift(MPI_Datatype & pop_type);  //PLUS LOADBALANCING HERE - blocking

	void PrintDensity(std::complex<double> trace);

	void Ms_RecvStat(MPI_Datatype & stat_type);  //get statistics - blocking (for master)
	void Ms_SendShift(MPI_Datatype & pop_type); //calculate and send stat. + later load balance here - non-blocking
	void Ms_CalcShift();
	void Ms_WriteStat(std::string filename, int time); //write the statistics into file
	void Ms_WriteStat_Support(std::string filename, int time); //write the statistics into file

	void Generate_Determinant(RepType det, int nbr_re, int nbr_im, System_Model & sys);
	void Generate_Connection_Col(Determinant_properties &properties, std::bitset<Bl> &det, std::bitset<Bl> &st, int & pos, int & ch_pos, System_Model & sys, std::complex<double> &hh, std::complex<double> &dhh);
	void Generate_Connection_Row(Determinant_properties &properties, std::bitset<Bl> &det, std::bitset<Bl> &st, int & pos, int & ch_pos, System_Model & sys, std::complex<double> &hh, std::complex<double> &dhh);
	void Generate_Field_Col(Determinant_properties &properties, std::bitset<Bl> &st, int & pos, System_Model & sys, std::complex<double> &ff, int col);

	int Death(std::complex<double> & diag_ampl, int parent, int child, int walker, int parent_sign);

	void Cleaner(); //to clean certain variables at the end of the time cycle

	void Ms_Write_Reset(int iter);  //write out parameters for restart, assuming the same input file
	void Sl_Write_Reset();  // write out determinants and walkers for restart
	void Sl_Write_Reset_Reference(); // write out the reference wave function

	void TestPerformance(Determinant d, int walk_nbr);
//	void Ms_Read_Reset(std::string filename,std::vector<std::bitset<Bl>> & refdets, std::vector<int> & walkers);

//	void Ms_Restart(MPI_Datatype & pop_c, MPI_Datatype & spawn_type, MPI_Datatype & bits);
//	void Sl_Restart(MPI_Datatype & pop_c, MPI_Datatype & spawn_type, MPI_Datatype & bits, System_Model & sys);
};
