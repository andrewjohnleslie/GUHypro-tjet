/*!
 *  \brief      Executable for the GP optimisation test case
 *  \author     Alessandro Mogavero
 *  \n
 *  \email      alessandro.mogavero@strath.ac.uk
 *  \version    1.0
 *  \copyright  Copyright 2016 Alessandro Mogavero
 */
 /* This file is part of HyPro.
 *
 * HyPro is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HyPro is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HyPro.  If not, see <http://www.gnu.org/licenses/>. */

#include "GP.h"

#ifdef MATLAB

void mexFunction( int nlhs, mxArray *plhs[],
		int nrhs, const mxArray *prhs[] )
{
	/* Check for proper number of arguments. */
	if(nlhs>1) {
		mexErrMsgIdAndTxt( "HYPRO:invalidNumOutputs",
				"Too many output arguments");
	}

	int gen = 20;
	int gp = -1;
	bool best = true;
	int rmId = -1;
	if(nrhs>0){
		if((!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) ||
				(mxGetN(prhs[0])!=1 || mxGetM(prhs[0])!=1)){
			mexErrMsgIdAndTxt( "HYPRO:invalidIO",
					"Amax must be a non-complex scalar");
		}

		gen = (int)(*mxGetPr(prhs[0]) + 0.5);
	}
	if(nrhs>1){
		if((!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])) ||
				(mxGetN(prhs[1])!=1 || mxGetM(prhs[1])!=1)){
			mexErrMsgIdAndTxt( "HYPRO:invalidIO",
					"Amax must be a non-complex scalar");
		}

		gp = (int)(*mxGetPr(prhs[1]) + 0.5);
		best = false;
	}
	if(nrhs>2){
		if((!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2])) ||
				(mxGetN(prhs[2])!=1 || mxGetM(prhs[2])!=1)){
			mexErrMsgIdAndTxt( "HYPRO:invalidIO",
					"Amax must be a non-complex scalar");
		}

		rmId = (int)(*mxGetPr(prhs[2]) + 0.5);
	}

	nonLinear<systemModule>::maxIter_ = 100;
	// We set up a new-handler, because we might need a lot of memory,
	// and we don't know it's there.

	// Init GP system.
	GPInit (1, -1);
	// Init HyPro
	thermoState::warning_ = false;
	MyPopulation::registerOptimClasses();

	// Declare the GP Variables, set defaults and read configuration
	// file.  The defaults will be overwritten by the configuration file
	// when read.  If it doesn't exist, the defaults will be written to
	// the file.
	GPConfiguration config (cout, "GP.ini", configArray);

	OptimSystem* sys;

	try{
		feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
		sys = inspect(gen, gp, true, rmId, 0, 0);
		fedisableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
	}catch (std::exception& e){
		fedisableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
		mexErrMsgIdAndTxt( "HYPRO:error", e.what());
	}

	plhs[0]  = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
	*reinterpret_cast<propModule**>(mxGetPr(plhs[0])) = (propModule*)sys;
}

#else

///Main test function
/**
 * SYNOPSIS
 * ========
 * \b GP [-out outLev] [-append] [-prevSeed] \n
 * \b GP -gen genID -gp best | gpID [-rm rmID]
 *
 * DESCRIPTION
 * ===========
 * The first syntax run a GP simulation.
 * The second syntax inspect preexisting data
 *
 * OPTIONS
 * =======
 * run options
 * -----------
 * <b>-out outLev</b> level of output verbosity. If 2 save details at each generation.
 * 	If 1 print advancement of optimisation. If 0 silent.\n
 * <b>-append</b> append new result to pre-existing generation data \n
 * <b>-prevSeed</b> use seed stored in seed.dat during a previous run
 *
 * inspect options
 * ---------------
 * <b>-gen genID</b> generation number of the population to inspect\n
 * <b>-gp best|gpID</b> select the GP individual to inspect.\n
 * 	If best select the best individual in the population.
 * 	If gpID number is specified select the corresponding GP individual.\n
 * <b>-rm rmId</b> ID of the module to be removed from the system module, for testing purposes.
 */
int main (int argc, char* argv[])
{
	inDir = "/home/mark/git/GUHypro/test/GP/resources/supersonic"; //Read inputs from work directory
//	inDir = "."; //Read inputs from work directory

	//Read standard input options
	int genID = -1;
	int gpID = -1;
	bool best = false;
	int rmId = -1;
	int outLev = 2;
	bool append = false;
	bool prevSeed = false;
//	std::cout<<"\nArgc is: " << argc;
//	std::cout<<"\nArgv is: " << argv[1] <<"\n";

	for(int i=1; i<argc; i+=2){
		if(strcmp(argv[i],"-gen")==0){
			sscanf(argv[i+1],"%d",&genID);
		}else if(strcmp(argv[i],"-gp")==0){
			if(strcmp(argv[i+1],"best")==0){
				best = true;
			}else{
				sscanf(argv[i+1],"%d",&gpID);
			}
		}else if(strcmp(argv[i],"-rm")==0){
			sscanf(argv[i+1],"%d",&rmId);
		}else if(strcmp(argv[i],"-out")==0){
			sscanf(argv[i+1],"%d",&outLev);
		}else if(strcmp(argv[i],"-append")==0){
			append = true;
		}else if(strcmp(argv[i],"-prevSeed")==0){
			prevSeed = true;
		}else{
			std::cout << (argv[i]) << std::endl;
			std::cerr << "Error: wrong input." << std::endl;
			exit(1);
		}
	}

	if(genID>-1){
		prevSeed = true; //Prevent changing the feed file in inspect mode
	}
	setup(prevSeed);

	if(genID!=-1 && (gpID!=-1 || best)){
		int Argc = 1;
		char* Argv2 = argv[0];
		char** Argv = &Argv2;
		inspect(genID, gpID, best, rmId, Argc, Argv);
	}else{
//		run(outLev,append);
	}

	return 0;
}

#endif
