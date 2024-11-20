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


int main(int argc, char *argv[]) {
	//Gas model construction
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

	if(argc<2){
		throw std::runtime_error("Error: 1 input requested.");
	}

	thermoState::warning_ = false;


	//Read HyPro model
	std::string fileName(argv[1]);
	std::ifstream fil(fileName);
	if(fil.fail()){
		throw std::runtime_error("Error: cannot open file " + fileName);
	}
	std::map<unsigned, std::shared_ptr<Node> > nodeMap;
	std::map<unsigned, std::shared_ptr<propModule> > moduleMap;
	systemModule::NodePtr N1 = Node::unserialize(fil,specieGas);
	systemModule::NodePtr N2 = Node::unserialize(fil,specieGas);
	nodeMap[N1->ID_] = N1;
	nodeMap[N2->ID_] = N2;

	systemModule::ModulePtr sys = propModule::unserialize(fil,nodeMap,moduleMap);
	fil.close();

	//Wedge Modification of FreeStream Conditions

	//Freestream conditions, Alt = 54,500ft

	double p = 9341.65;
	double T = 216.65;
	double Mach = 3.5;

	Node freeStream(specieGas);
	freeStream.p_ = p;
	freeStream.T_ = T;
	List<double> X(4,0.0);
	X[1] = 0.209;   //Air
	X[2] = 1-0.209;
	freeStream.X(X);
	freeStream.U_ = Mach*freeStream.a();
	freeStream.Name_ = "Free Stream";
	Node freeStreamMod(specieGas);
	freeStreamMod.Name_ = "Pre Intake";

	Wedge forebody("forebody",&freeStream,&N1);
	forebody.delta_ = 8.0;
	forebody.calculate();

	sys->unreduce();
	sys->calculate();
	systemModule* sysC = (systemModule*)sys.get();
	sysC->printOut(std::cout);

	std::cout << "T0:" << std::endl;

	for(std::size_t i=0; i<sysC->nodes_.size(); i++){
	double H = sysC->nodes_[i]->H();
	double T0 = sysC->nodes_[i]->gas().TH(H,300);
	std::cout << sysC->nodes_[i]->Name_ << " " << T0 << std::endl;
		}

	double H = sysC->N2().H();
	double T0 = sysC->N2().gas().TH(H,300);
	std::cout << "Nozzle Exit " << T0 << std::endl << std::endl;


	std::cout << "Thrust = " << sysC->thrust() << std::endl;

	//Save model
	std::ofstream fil1(fileName + ".out");
	sys->serialize(fil1);
	fil1.close();

	return 1;
}
