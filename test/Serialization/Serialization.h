/*!
 *  \brief      Verification of the serialization procedure
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
#include <time.h>
#include <map>

#include "core/Node.h"

#include <IOdictionary.H>
#include <IOobject.H>
#include <Time.H>
#include <argList.H>
#include <specieThermo.H>
#include <janafThermo.H>
#include <perfectGas.H>

using namespace Foam;

#include <mixture.h>
#include "solvers/isentropicDuct.h"
#include "core/systemModule.h"
#include "combustion/EffPhiComb.h"
#include "injectors/InjectionPhi.h"
#include <AdaptedThroatInlet.h>
#include "Wedge.h"
#include "Rocket.h"
#include "Mixer.h"
#include <AdaptedNozzle.h>
#include "core/MultiModeMachP0.h"
#include <AdaptedInlet.h>
#include "NeutralLink.h"
#include <chokedConv.h>
#include <solvers/Friction.h>
#include <solvers/EffMachIsenDuct.h>
#include <injectors/InjectionPlate.h>
#include <AdaptedThroatInlet.h>

std::unique_ptr<systemModule> createModel(){

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

	mixture< specieThermo< janafThermo<perfectGas> > >& specieGas =
			*new mixture< specieThermo< janafThermo<perfectGas> > >(a,species);

	std::vector<double> X(4,0.0);
	X[1] = 0.209;   //Air
	X[2] = 1-0.209;

	thermoState::warning_ = false;

	Node& cfdInlet(*new Node(specieGas));
	cfdInlet.Name_ = "CFD Inlet";
	cfdInlet.p_ = 1e5;
	cfdInlet.T_ = 300.0;
	cfdInlet.X(X);
	Node& exit(*new Node(cfdInlet));
	exit.Name_ = "Nozz exit";

	std::unique_ptr<systemModule> sysPtr(new systemModule("SERJ",&cfdInlet,&exit));
	systemModule& system(*sysPtr);
	//system.verbosity_ = 0;

	const int numNodes = 10;
	double A[numNodes] =   {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};

	cfdInlet.A(A[0]);
	system.initNodes(numNodes,cfdInlet);
	for(int i=0; i<numNodes; i++){
	    system.nodes_[i]->A(A[i]);
	}
	exit.A(3.5);

	// Pre-primary
	systemModule::ModulePtr mod1(new EffIsenDuct("mod1",&system.N1(),system.nodes_[0].get()));
	mod1->N2().Name_ = "mod1-out";
	system.add(mod1);

	systemModule::ModulePtr mod2(new Mixer("mod2",system.nodes_[0].get(),
			system.nodes_[1].get(),system.nodes_[2].get()));
	mod2->N2().Name_ = "mod2-out";
	Mixer* primaryC = (Mixer*)mod2.get();
	primaryC->teta_ = 0;
	primaryC->eff_ = 0.7;
	system.add(mod2);

	// Nozzle
	systemModule::ModulePtr mod3(new chokedConv("mod3",system.nodes_[2].get(),
			system.nodes_[3].get(),chokedConv::CALCULATED));
	system.add(mod3);
	mod3->N2().Name_ = "mod3-out";

	systemModule::ModulePtr mod4(new Friction("mod4",system.nodes_[3].get(),
				system.nodes_[4].get()));
	system.add(mod4);
	mod4->N2().Name_ = "mod4-out";
	Friction* mod4C = (Friction*)mod4.get();
	mod4C->Cf_ = 0.6;
	mod4C->D_ = 2.0;
	mod4C->L_ = 1.6;

	systemModule::ModulePtr mod5(new EffPhiComb<Friction>("mod4",system.nodes_[4].get(),
					system.nodes_[5].get()));
	system.add(mod5);
	mod5->N2().Name_ = "mod5-out";
	EffPhiComb<Friction>* mod5C = (EffPhiComb<Friction>*)mod5.get();
	mod5C->Pcoeff_[0] = 3.0;
	mod5C->Rcoeff_[2] = 2.4;

	systemModule::ModulePtr mod6(new EffMachIsenDuct("mod6",system.nodes_[5].get(),system.nodes_[6].get()));
	mod6->N2().Name_ = "mod6-out";
	system.add(mod6);

	systemModule::ModulePtr mod7(new InjectionPhi<Friction>("mod7",system.nodes_[6].get(),system.nodes_[7].get()));
	mod7->N2().Name_ = "mod7-out";
	system.add(mod7);

	systemModule::ModulePtr mod8(new InjectionPlate("mod8",system.nodes_[7].get(),system.nodes_[8].get()));
	mod8->N2().Name_ = "mod8-out";
	system.add(mod8);

	systemModule::ModulePtr mod9(new AdaptedThroatInlet("mod8",system.nodes_[8].get(),system.nodes_[9].get()));
	mod9->N2().Name_ = "mod9-out";
	system.add(mod9);

	return sysPtr;
}

void deleteAll(systemModule& sys){
	const mixture< specieThermo< janafThermo<perfectGas> > >* gas = &sys.N1().specieGas();
	delete(&sys.N1());
	delete(&sys.N2());
	delete(gas);
}
