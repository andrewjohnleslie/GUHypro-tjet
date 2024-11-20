/*!
 *  \brief      Validation test case based on the Hyperion propulsion system
 *  \details    Hyperion is a space plane vehicle developed at the GeorgiaTech.
 *              Details of the validation are available in the
 *              <a href="http://oleg.lib.strath.ac.uk/R/?func=dbin-jump-full&object_id=26892">PhD thesis of Alessandro Mogavero</a>
 *
 *              @param[in] p free stream static pressure
 *              @param[in] T free stream static temperature
 *              @param[in] Mach free stream Mach number
 *              @param[in] mode propulsion mode number, -1 if auto mode selection is wanted
 *              @param[in] verbosity 0 if no output is wanted, 1 only for node output, 2 for additional output during calculation
 *              @param[out] Thrust net thrust
 *              @param[out] Mfr mass flow rate of propellant for each chemical species
 *              @param[out] out detailed output at each engine station
 *              @return a pointer to the engine system model
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

#include "Node.h"

// #include <IOdictionary.H>
// #include <IOobject.H>
// #include <Time.H>
// #include <argList.H>
// #include <specieThermo.H>
// #include <janafThermo.H>
// #include <perfectGas.H>
// using namespace Foam;

// #include <mixture.h>
#include "isentropicDuct.h"
#include "systemModule.h"
#include "EffComb.h"
#include "InjectionPhi.h"
#include <AdaptedThroatInlet.h>
#include "Wedge.h"
#include "Rocket.h"
#include "Mixer.h"
#include <AdaptedNozzle.h>
#include "MultiModeMachP0.h"
#include <AdaptedInlet.h>
#include "NeutralLink.h"
#include "EffIsenDuct.h"
#include "Friction.h"

using namespace hypro;
systemModule* Hyperion(const double& p, const double& T, const double& Mach, const int mode,
		const int verbosity, double& Thrust, std::vector<double>& Mfr, Collection::OutType& out){

	double foot =  0.3048; //conversion feet to meters
	double lbm = 0.45359; //conversion libre to kg

	// // Thermal model
	// Time runTime(QUOTE(HyProRescource),"");

	// IOdictionary a
	// (
	// 		IOobject
	// 		(
	// 				"thermophysicalProperties",
	// 				"",
	// 				runTime,
	// 				IOobject::MUST_READ,
	// 				IOobject::NO_WRITE,
	// 				false
	// 		)
	// );
	// wordList species(4,"H2");
	// species[1] = "O2";
	// species[2] = "N2";
	// species[3] = "H2O";

	// mixture< specieThermo< janafThermo<perfectGas> > > specieGas(a,species);

	// Node freeStream(specieGas);
	// freeStream.p_ = p;
	// freeStream.T_ = T;
	// List<double> X(4,0.0);
	// X[1] = 0.209;   //Air
	// X[2] = 1-0.209;
	// freeStream.X(X);
	// freeStream.U_ = Mach*freeStream.a();
	// freeStream.Name_ = "Free Stream";

//	Node freeStream("/usr/share/cantera/data/gri30_highT.yaml", "gri30");
	Node freeStream("/home/robert/repos/HyPro/test/Hyperion/gri30_highT2.xml", "gri30_simple");
    double pStream = p;
	if (p < 1){
		pStream = 1;
		std::cout << "P too low" << std::endl;
	}

	Cantera::compositionMap initialComposition;
	initialComposition.emplace("O2",0.209);
	initialComposition.emplace("N2",1-0.209);
	
	freeStream.setTPX(T,pStream,initialComposition);
	vector<double> empty_thermos(freeStream.X().size());
	freeStream.setU(Mach*freeStream.geta());
	freeStream.Name_ = "Free Stream";

	thermoState::warning_ = false;

	// Node freeStreamMod(specieGas);
	// freeStreamMod.Name_ = "Pre Intake";
	// Node exit(specieGas);
	// exit.Name_ = "Nozzle Exit";

	Node& freeStreamMod = *new Node("/home/robert/repos/HyPro/test/Hyperion/gri30_highT2.xml", "gri30_simple");
	freeStreamMod.Name_ = "Pre Intake";
	Node& exit = *new Node("/home/robert/repos/HyPro/test/Hyperion/gri30_highT2.xml", "gri30_simple");
	exit.Name_ = "Nozzle Exit";

	Wedge forebody("forebody",&freeStream,&freeStreamMod);
	forebody.delta_ = 8.0;
	forebody.calculate();

	MultiModeMachP0 multiSystem("system",&freeStreamMod,&exit);
	multiSystem.verbosity_ = verbosity - 1;
	systemModule& system(multiSystem.add("Ejector Mode"));
	system.verbosity_ = verbosity - 1;

	//Nodes                 N1_    I0    I1       0     1       2     3      4        5      6     7     N2_
	const int numNodes = 12;
	double A[numNodes] =   {27.0, 27.0, 0.25*27, 8.24, 11.25, 11.25, 0.0, 11.25-8.24, 22.5, 22.5, 22.5, 44.3};
	/*name = {'Pre Intake','Intake','Throat','Pinch Point','Primary','Mixer End',...
	    'Dummy','Rocket Outlet','End Diffuser','Injection','End Chamber','Nozzle'};*/
	double foot2 = std::pow(foot,2.0);
	freeStreamMod.A(A[0]*foot2);
	system.initNodes(8,freeStreamMod);
	for(int i=0; i<8; i++){
	    system.nodes_[i]->A(A[i+3]*foot2);
	}
	exit.A(A[numNodes-1]*foot2);
	exit.Amax_ = A[numNodes-1]*foot2;

	//Inlet
	systemModule::ModulePtr inlet(new AdaptedThroatInlet("Inlet",&system.N1(),system.nodes_[0].get()));
	inlet->N2().Name_ = "Pinch Point";
	AdaptedThroatInlet* inletC = (AdaptedThroatInlet*)inlet.get();
	inletC->nodes_[0]->A(A[1]*foot2);
	inletC->nodes_[0]->Amax_ = inletC->nodes_[0]->A();
	inletC->nodes_[0]->Name_ = "Intake";
	inletC->nodes_[1]->A(A[2]*foot2);
	inletC->nodes_[1]->Amax_ = inletC->nodes_[1]->A();
	inletC->nodes_[1]->Name_ = "Throat";
	inlet->verbosity_ =0;
	inlet->chokingFeedback_ = inlet;
	system.add(inlet);

	// Primary - ejector
	systemModule::ModulePtr postPinch(new isentropicDuct(
			"Post pinch point",&inlet->N2(),system.nodes_[1].get()));
	system.add(postPinch);

	systemModule::ModulePtr rocket(new Rocket("Rocket",system.nodes_[3].get(),system.nodes_[4].get()));
	rocket->N2().Name_ = "Rocket exit";
	Rocket* rocketC = (Rocket*)rocket.get();
	// List<double> Xr(4,0.0);
	std::vector<double> Xr = empty_thermos;
	// Xr[3] = 1;
	Xr[2] = 1;
	rocketC->N2().X(Xr);
	rocketC->Me_ = 6.0;
	rocketC->mfrMax_ = 216*lbm;
	rocketC->unreduce();
	rocketC->Isp_ = 462;
	rocketC->unreduce();
	system.add(rocket);

	systemModule::ModulePtr primary(new Mixer("primary",
			&postPinch->N2(),system.nodes_[2].get(),&rocket->N2()));
	primary->N1().Name_ = "Primary In";
	primary->N2().Name_ = "Primary Out";
	Mixer* primaryC = (Mixer*)primary.get();
	primaryC->verbosity_ = 0;
	primaryC->chokingFeedback_ = inlet;
	primaryC->eff_ = 0.9;
	system.add(primary);

	// Secondary - Post combustor
	systemModule::ModulePtr diffuser(new EffIsenDuct("diffuser",system.nodes_[2].get(),system.nodes_[5].get()));
	diffuser->N2().Name_ = "Diffuser End";
	EffIsenDuct* diffuserC = (EffIsenDuct*)diffuser.get();
	diffuserC->etap_ = 0.92;
	diffuserC->etaT_ = 1.0;
	system.add(diffuser);

	systemModule::ModulePtr secondary(new InjectionPhi<balanceMach>("secondary",system.nodes_[5].get(),system.nodes_[6].get()));
	secondary->N2().Name_ = "Injection";
	InjectionPhi<balanceMach>* secondaryC = (InjectionPhi<balanceMach>*)secondary.get();
	secondaryC->PhiMax_ = 1;
	secondaryC->Phi_ = 1;
	// List<double> Xf(4,0.0);
	std::vector<double> Xf = empty_thermos;
	Xf[0] = 1;
	secondaryC->Xf_ = Xf;
	secondaryC->Tf_ = 20;
	secondaryC->Rcoeff_[0] = 2;
	// secondaryC->Rcoeff_[1] = 1;
	secondaryC->Rcoeff_[1] = 1;
	secondary->chokingFeedback_ = secondary;
	system.add(secondary);

	systemModule::ModulePtr combustor(new EffComb<Friction>("combustor",system.nodes_[6].get(),system.nodes_[7].get()));
	combustor->N2().Name_ = "Chamber End";
	combustor->chokingFeedback_ = secondary;
	EffComb<Friction>* combustorC = (EffComb<Friction>*)combustor.get();
	combustorC->Rcoeff_[0] = 2;
	// combustorC->Rcoeff_[1] = 1;
	// combustorC->Pcoeff_[3] = 2;
	combustorC->Rcoeff_[1] = 1;
	combustorC->Pcoeff_[2] = 2;
	combustorC->eta_ = 0.8;
	combustorC->Cf_ = 0.01;
	combustorC->D_ = 1.631;
	combustorC->L_ = 1.0;
	system.add(combustor);

	// Nozzle
	systemModule::ModulePtr nozzle(new AdaptedNozzle<EffIsenDuct>("nozzle",system.nodes_[7].get(),&system.N2()));
	AdaptedNozzle<EffIsenDuct>* nozzleC = (AdaptedNozzle<EffIsenDuct>*)nozzle.get();
	nozzleC->choked_ = true;
	/* in case of inlet spilling (especially in supersonic) better to adapt the
		nozzle to the inlet conditions.*/
	nozzleC->freeStream_ = &inlet->N1();
	nozzleC->etap_ = 1.0;
	system.add(nozzle);

	// Ramjet mode
	systemModule& ramjet(multiSystem.add("Ramjet Mode",0));
	const systemModule::ModulePtr inlet1(ramjet.exchange<AdaptedThroatInlet>(inlet,"Inlet"));
	AdaptedThroatInlet* inlet1C = (AdaptedThroatInlet*)inlet1.get();
	inlet1C->nodes_[0]->A(A[1]*foot2);
	inlet1C->nodes_[0]->Amax_ = inlet1C->nodes_[0]->A();
	inlet1C->nodes_[0]->Name_ = "Intake";
	inlet1C->nodes_[1]->A(1.055*A[2]*foot2);
	inlet1C->nodes_[1]->Amax_ = inlet1C->nodes_[1]->A();
	inlet1C->nodes_[1]->Name_ = "Throat";
	inlet->verbosity_ =0;
	inlet->chokingFeedback_ = inlet;

	const systemModule::ModulePtr primary1(ramjet.exchange<EffIsenDuct>(primary,"primary"));
	primary1->chokingFeedback_ = inlet1;
	EffIsenDuct* primary1C = (EffIsenDuct*)primary1.get();
	primary1C->etap_ = 0.725;

	ramjet.remove(rocket);

	//Pure Rocket mode
	systemModule& rocketMode(multiSystem.add("Rocket Mode",0));
	const systemModule::ModulePtr primary2(rocketMode.exchange<isentropicDuct>(primary,"primary"));
	primary2->N1(rocket->N2());
	const systemModule::ModulePtr secondary2(rocketMode.exchange<NeutralLink>(secondary,"secondary"));
	secondary2->N2(combustor->N2());
	rocketMode.remove(inlet);
	rocketMode.remove(postPinch);
	rocketMode.remove(combustor);
	Node rocketDummyN1(rocketMode.N1());
	rocketMode.N1(rocketDummyN1);
	rocketDummyN1.A(0.0);

	multiSystem.rangesM_[0][0] = 0.0;
	multiSystem.rangesM_[0][1] = 3.0;
	multiSystem.rangesM_[1][0] = 3.0;
	multiSystem.rangesM_[1][1] = 7.0;
	multiSystem.rangesP0_[0][0] = 1e4;
	multiSystem.rangesP0_[1][0] = 5e5;

	if(mode==-1){
		multiSystem.calculate();
		Thrust = multiSystem.thrust();
		Mfr = multiSystem.propellantMfr();
	}else if(mode>=0){
		systemModule& selMode(multiSystem.modes_[mode]);
		selMode.verbosity_ = verbosity - 1;
		selMode.calculate();
		Thrust = selMode.thrust();
		Mfr = selMode.propellantMfr();

		Collection::OutType out1(selMode.out());
		for(std::size_t i=0; i<out1.size(); i++){
			out.push_back(out1[i]);
		}

		if(verbosity>0){
			selMode.printOut(std::cout);
		}
	}

	return new systemModule(system);
}
