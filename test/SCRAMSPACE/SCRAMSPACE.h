/*!
 *  \brief      Validation against CFD data of a scramjet test bed
 *  \details    Details of the validation are available in <a href="http://oleg.lib.strath.ac.uk/R/?func=dbin-jump-full&object_id=26892">
 *              "Toward automated design of Combined Cycle Propulsion." by Alessandro Mogavero</a>
 *
 *              @param[in] p inlet static pressure
 *              @param[in] T inlet static temperature
 *              @param[in] U inlet speed
 *              @param[in] combEff combustion efficiency
 *              @param[in] verbosity 0 if no output is wanted, 1 only for node output, 2 for additional output during calculation
 *              @param[out] Thrust net thrust
 *              @param[out] Mfr mass flow rate of propellant for each chemical species
 *              @param[out] out detailed output at each station of the engine
 *              @return pointer to the system module
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

#pragma once

#define DO_QUOTE(X)       #X
#define QUOTE(X)          DO_QUOTE(X)

#include <iostream>
#include <fstream>
#include <time.h>

//#include <IOdictionary.H>
//#include <IOobject.H>
//#include <Time.H>
//#include <argList.H>
//#include <specieThermo.H>
//#include <janafThermo.H>
//#include <perfectGas.H>
//using namespace Foam;
//#include <mixture.h>

#include "core/Node.h"
#include "core/systemModule.h"
#include "solvers/isentropicDuct.h"
#include "solvers/Friction.h"
#include "Fanno.h"
#include "thermo/thermoState.h"

#include "solvers/EffIsenDuct.h"
#include "combustion/EffComb.h"
#include "injectors/injection.h"
#include "InjectionPhi.h"

#include "cantera/thermo/ThermoPhase.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/Phase.h"

using namespace hypro;

systemModule* SCRAMSPACE(const double& p, const double& T, const double& U, const double& combEff, const int& verbosity,
		double& Thrust, std::vector<double>& Mfr, Collection::OutType& out){

	Cantera::compositionMap initialComposition;
	initialComposition.emplace("O2",0.209);
	initialComposition.emplace("N2",1-0.209);

	//Node& freeStream = *new Node("/usr/share/cantera/data/gri30_highT.xml", "gri30");
	Node& freeStream = *new Node("/usr/share/cantera/data/gri30_highT.yaml", "gri30");
	freeStream.setTPX(T,p,initialComposition);
    freeStream.setU(U);
    freeStream.Name_ = "Free Stream";

    vector<double> empty_thermos = freeStream.X();
    empty_thermos[3] = 0;
    empty_thermos[47] = 0;

	thermoState::warning_ = false;

	Node& exit = *new Node(freeStream);
	exit.Name_ = "Nozzle Exit";

	systemModule system("system",&freeStream,&exit);
	system.verbosity_ = verbosity - 1;

	//Nodes                    N1_         0       1       2         3         N2_
	const int numNodes = 6;
	double A[numNodes] =   {1*0.04385, 1*0.013, 1*0.013, 1*0.013, 1*0.013, 1*0.04385};

	freeStream.A(A[0]);
	system.initNodes(4,freeStream);
	for(int i=0; i<4; i++){
	    system.nodes_[i]->A(A[i+1]);
	}
	exit.A(A[numNodes-1]);
	exit.Amax_ = A[numNodes-1];

	//Inlet
	systemModule::ModulePtr inlet(new EffIsenDuct("Inlet",&system.N1(),system.nodes_[0].get()));
	inlet->N2().Name_ = "Inlet subsonic";
	EffIsenDuct* inletC = (EffIsenDuct*)inlet.get();
	inletC->etap_ = 0.45;
	inletC->choked_ = false;
	system.add(inlet);

//	//Injection
//	systemModule::ModulePtr inj(new injection(
//				"Injection",system.nodes_[0].get(),system.nodes_[1].get()));
//	inj->N2().Name_ = "Injection end";
//	injection* injC = (injection*)inj.get();
//    std::vector<double> Xf = empty_thermos;
//    Xf[0] = 0.25;
//    injC->Xf_ = Xf;
//    system.add(inj);

	//Injection
	systemModule::ModulePtr inj(new InjectionPhi<balanceMach>(
				"Injection",system.nodes_[0].get(),system.nodes_[1].get()));
	inj->N2().Name_ = "Injection end";
	InjectionPhi<balanceMach>* injC = (InjectionPhi<balanceMach>*)inj.get();
    std::vector<double> Xf = empty_thermos;
    injC->Xf_ = Xf;
    injC->Xf_[0] = 1;
    injC->Phi_ = 1.0;
	injC->PhiMax_ = 1.0;
	injC->Tf_ = 350;
    injC->Rcoeff_[0] = 2;
    injC->Rcoeff_[3] = 1;
    system.add(inj);

	// isolator
	systemModule::ModulePtr isolator(new Friction(
			"Isolator",system.nodes_[1].get(),system.nodes_[2].get()));
	isolator->N2().Name_ = "Isolator end";
	Friction* isolatorC = (Friction*)isolator.get();
	isolatorC->L_ = 0.15;
	isolatorC->D_ = 4*0.026*1/(2*1);
	isolatorC->Cf_ = 0.0035;
	system.add(isolator);

	systemModule::ModulePtr combustor(new EffComb<Friction>("combustor",system.nodes_[2].get(),system.nodes_[3].get()));
	combustor->N2().Name_ = "Chamber End";
	EffComb<Friction>* combustorC = (EffComb<Friction>*)combustor.get();
	combustorC->Rcoeff_[0] = 2;
	combustorC->Rcoeff_[3] = 1;
	combustorC->Pcoeff_[5] = 2;
	combustorC->eta_ = combEff;
	combustorC->Cf_ = 0.0035;
	combustorC->D_ = 4*0.026*1/(2*1);
	combustorC->L_ = 0.25;
	system.add(combustor);

	// Nozzle
	systemModule::ModulePtr nozzle(new isentropicDuct("nozzle",system.nodes_[3].get(),&system.N2()));
	system.add(nozzle);

	system.calculate();
    Thrust = system.thrust();
	Mfr = system.propellantMfr();
	out = system.out();
	if(verbosity>0){
		system.printOut(std::cout);
	}
//    std::cout << system.nodes_[1].get()->report() << std::endl;
//	std::cout << system.nodes_[3].get()->report() << std::endl;
//	std::cout << system.N2().report() << std::endl;
	return new systemModule(system);
}
