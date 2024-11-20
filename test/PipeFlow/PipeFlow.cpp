/*!
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

#define WM_DP
#define NoRepository

#define DO_QUOTE(X)       #X
#define QUOTE(X)          DO_QUOTE(X)

#include <iostream>
#include <fstream>
#include <time.h>

#include "core/Node.h"

#include <IOdictionary.H>
#include <IOobject.H>
#include <Time.H>
#include <argList.H>
#include <specieThermo.H>
#include <janafThermo.H>
#include <perfectGas.H>
#include <perfectGas.C>
using namespace Foam;

#include <mixture.h>
#include "solvers/isentropicDuct.h"
#include "core/systemModule.h"
#include "combustion/Combustion.h"
#include "injectors/InjectionPhi.h"
#include <AdaptedThroatInlet.h>
#include "Wedge.h"
#include "Rocket.h"
#include "Mixer.h"
#include <AdaptedNozzle.h>
#include "core/MultiModeMachP0.h"
#include <AdaptedInlet.h>
#include "NeutralLink.h"
#include "solvers/EffIsenDuct.h"
#include "solvers/Friction.h"
//#include "ISAModel.h"
#include "Fanno.h"


#ifdef MATLAB
#include "mex.h"
#include <fenv.h>
#endif

systemModule* Pipeflow(const double& p, const double& T, const double& Mach, const int mode,
		double& Thrust, Foam::List<double>& Mfr, Collection::OutType& out){

	// Thermal model
	Time runTime(QUOTE(HyProRescource),"");

	IOdictionary a
	(
			IOobject
			(
					"thermophysicalProperties",
					"",
					runTime,
					IOobject::MUST_READ,
					IOobject::NO_WRITE,
					false
			)
	);
	wordList species(4,"H2");
	species[1] = "O2";
	species[2] = "N2";
	species[3] = "H2O";

	mixture< specieThermo< janafThermo<perfectGas> > > specieGas(a,species);

	Node freeStream(specieGas);
	freeStream.p_ = p;
	freeStream.T_ = T;
	List<double> X(4,0.0);
	X[1] = 0.209;   //Air
	X[2] = 1-0.209;
	freeStream.X(X);
	freeStream.U_ = Mach*freeStream.a();
	freeStream.Name_ = "Pipe Inlet";

	thermoState::warning_ = false;

	Node exit(specieGas);
	exit.Name_ = "Pipe Exit";


	MultiModeMachP0 multiSystem("system",&freeStream,&exit);
	multiSystem.verbosity_ = 1;
	systemModule& system(multiSystem.add("Fanno Mode"));
	system.verbosity_ = 0;

	//Nodes                  N1_   N2_
	const int numNodes = 2;
	double A[numNodes] =   {1.0,  1.0};
	freeStream.A(A[0]);
	exit.A(A[1]);

	propModule::chokingMach_ = 1.0;


	// Select Friction class to find normal Hypro results, or Fanno to solve analytically
		// Both classes require definition of Cf, D and L and this rematins untouched when choosing which class to use, only swap the first 4 lines

	//systemModule::ModulePtr fanno(new Fanno("fanno",&system.N1(),&system.N2()));
	systemModule::ModulePtr fanno(new Friction("fanno",&system.N1(),&system.N2()));
	//Fanno* fannoC = (Fanno*)fanno.get();
	Friction* fannoC = (Friction*)fanno.get();

	fanno->N2().Name_ = "Fanno Pipe End";
	fannoC->Cf_ = 0.01;
	fannoC->D_ = 1.128; 				// D here was caluclated by hand, but relates to an area of 1.0 sq m
	//fannoC->L_ = 6;
	system.add(fanno);


	//Iterative Length Solver
	//Used when trying to solve for an unknown length (searching for the length and physical conditions with a choked exit)
	//To Use iterative solver to find length that will choke flow at exit, comment out line "fannoC->L_ = double;" above



	double toll_ = 1e-6;
	double maxIter_ = 1000;


	double x1 = 0.0;
	double x2 = 25;



	fannoC->L_ = x1;
	double y1 = fanno->calculate();

	fannoC->L_ = x2;
	double y2 = fanno->calculate();

	while(y1*y2>0){
			//throw std::runtime_error("Error: function evaluations at starting points must be of opposite signs.");
			x1 = x2;
			fannoC->L_ = x1;
				y1 = fanno->calculate();

			x2 = 1.5*x2;
			fannoC->L_ = x2;
				y2 = fanno->calculate();
		}

		double r = 2*toll_;
		int i = 0;
		double x;
		while(r>toll_ and i<maxIter_){
			x = 0.5*(x1 + x2);

			fannoC->L_ = x;
			double y = fanno->calculate();


			if(y*y1<0){
				y2 = y;
				x2 = x;
			}
			else{
				y1 = y;
				x1 = x;
			}

			r = std::abs(x1-x2);
			i++;
		}

		if(i==maxIter_){
			std::cout << "Warning: maximum number of iteration reached." << std::endl;
		}

	fannoC->L_=x-0.000005;


	// End of Iterative Solver


	if(mode==-1){
		multiSystem.calculate();
		Thrust = multiSystem.thrust();
		Mfr = multiSystem.propellantMfr();
	}else{
		systemModule& selMode(multiSystem.modes_[mode]);
		selMode.verbosity_ = 0;
		selMode.calculate();
		Thrust = selMode.thrust();
		Mfr = selMode.propellantMfr();
	 	selMode.printOut(std::cout);

		std::vector<double> outTmp;
		outTmp.reserve(10);
		outTmp.push_back(freeStream.A());
		outTmp.push_back(freeStream.Amax_);
		outTmp.push_back(freeStream.Amin_);
		outTmp.push_back(freeStream.p_);
		outTmp.push_back(freeStream.T_);
		outTmp.push_back(freeStream.U_);
		outTmp.push_back(freeStream.X()[0]);
		outTmp.push_back(freeStream.X()[1]);
		outTmp.push_back(freeStream.X()[2]);
		outTmp.push_back(freeStream.X()[3]);
		out.push_back(std::make_pair(freeStream.Name_,outTmp));


		/*double s1 = 0.5 * system.N1().rho() * std::pow(system.N1().U_,2);
		double s2 = 0.5 * system.N2().rho() * std::pow(system.N2().U_,2);

		double friction = 0.01 * ((s1 + s2)/2) * (Length * 3.14 * 1.128);
		double difference = system.N2().I() - system.N1().I();

		std::cout <<  difference  <<  "    " << friction <<  "      " << friction/difference << std::endl ;*/
		//std::cout << s2 << std::endl;

		//double s1 = system.N1().gamma();
		//double s2 = system.N2().rho();

		//std::cout << "Exit Density:  " << s2 << std::endl;
		//std::cout << "Inlet Stagnation Temp: " << system.nodes_[0]->T0() << std::endl;
		//std::cout << "Exit Stagnation Temp: " << system.N2().T0() << std::endl;

		std::cout << x << std::endl;

		Collection::OutType out1(selMode.out());
		for(std::size_t i=0; i<out1.size(); i++){
			out.push_back(out1[i]);
		}
	}


	return new systemModule(system);
}




int main(int argc, char *argv[]) {

	/*int Machiterations = 1;
	double MachRange[1] = {2};*/

	int Machiterations = 21;
	double MachRange[21] = {3.0, 2.8, 2.6, 2.4, 2.2, 2.0, 1.8, 1.6, 1.4, 1.2, 1.1, 0.999, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1};





 		double p = 101320;
	 	double T = 1000;

	 	//std::cout  << std::endl << std::endl <<  "Pressure: " << p <<" Temperature: " << T << std::endl << "Local Speed of Sound: " << std::pow(1.4*287*T,0.5);



	 	for(int i=0; i<Machiterations; i++){

	 		double Mach = MachRange[i];

	 		std::cout  << std::endl << "NEW" << std::endl <<  "Mach Number: " << Mach << std::endl << std::endl;



	 			 		double Thrust;
	 			 		Foam::List<double> Mfr;
	 			 		Collection::OutType out;

	 			 		Pipeflow(p,T,Mach,0,Thrust,Mfr,out);




		}


 	return 1;
}
