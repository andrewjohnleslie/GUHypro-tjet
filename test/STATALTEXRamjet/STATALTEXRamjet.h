/*!
 *  \brief      Validation against ramjet flight data
 *  \details    For more details see section 4.6 in <a href="http://oleg.lib.strath.ac.uk/R/?func=dbin-jump-full&object_id=26892">
 *              "Toward automated design of Combined Cycle Propulsion." by Alessandro Mogavero</a>
 *
 *              @param[in] p free stream static pressure
 *              @param[in] T free stream static pressure
 *              @param[in] Mach free stream Mach number
 *              @param[in] verbosity control the printed output while running
 *              @param[in] STX flight number
 *              @param[out] Thrust net thrust
 *              @param[out] Mfr mass flow rate of propellants for each chemical species
 *              @param[out] out detailed output for each engine section
 *              @return[out] pointer to the system module
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

#include "core/Node.h"

//#include <IOdictionary.H>
//#include <IOobject.H>
//#include <Time.H>
//#include <argList.H>
//#include <specieThermo.H>
//#include <janafThermo.H>
//#include <perfectGas.H>
//using namespace Foam;

//#include <mixture.h>
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
#include "solvers/EffIsenDuct.h"
#include "solvers/Friction.h"
#include <FixedInlet.h>
#include <ZitaInlet.h>

#include "cantera/thermo/ThermoPhase.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/Phase.h"

hypro::systemModule* Ramjet(const double& p, const double& T, const double& Mach, const int verbosity, const int STX,
		double& Thrust, std::vector<double>& Mfr, hypro::Collection::OutType& out){

	double foot =  1; //conversion feet to meters

	/*// Thermal model
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

	wordList species(5,"JETA");
	species[1] = "O2";
	species[2] = "N2";
	species[3] = "H2O";
	species[4] = "CO2";


	mixture< specieThermo< janafThermo<perfectGas> > > specieGas(a,species);
*/

    std::string directory(QUOTE(STXRescource));

    Cantera::compositionMap initialComposition;
    initialComposition.emplace("O2",0.209);
    initialComposition.emplace("N2",1-0.209);

    hypro::Node& freeStream = *new hypro::Node(directory + "/gri30_highT_modified_jetA.yaml", "gri30");
    freeStream.setTPX(T,p,initialComposition);
    freeStream.setU(Mach*freeStream.geta());
	freeStream.Name_ = "Free Stream";

    vector<double> empty_thermos = freeStream.X();
    empty_thermos[3] = 0;
    empty_thermos[47] = 0;

	hypro::thermoState::warning_ = false;

	hypro::Node& exit = *new hypro::Node(freeStream);
	exit.Name_ = "Nozzle Exit";

	hypro::MultiModeMachP0 multiSystem("system",&freeStream,&exit);
	multiSystem.verbosity_ = 1;
	hypro::systemModule& system(multiSystem.add("Ramjet"));
	system.verbosity_ = 0;

	//Nodes                  N1_      I0      0       1       2       3       4      5     N2_
	const int numNodes = 9;
	double A[numNodes] = {0.018,    0.0180, 0.0379, 0.058,  0.058, 0.0934, 0.0934,  0.0,   0.0};
	double fuelMfrMax;
	if(STX==5){
		A[7] = 0.0415;
		A[8] = 0.1308;
		fuelMfrMax = 1.1;
	}else if(STX==6){
		A[7] = 0.0394;
		A[8] = 0.1576;
		fuelMfrMax = 1.2;
	}else if(STX==9){
		A[7] = 0.0394;
		A[8] = 0.1576;
		fuelMfrMax = 1.3;
	}else{
		throw std::logic_error("Error: STX not recognized");
	}

	/*name = {'Pre Intake','After Intake', 'Injection Outlet', 'Pre-Combustor' , 'Nozzle Inlet', 'Nozzle Exit'};*/
	double foot2 = std::pow(foot,2.0);
	freeStream.A(A[0]*foot2);
	system.initNodes(6,freeStream);
	for(int i=0; i<6; i++){
	    system.nodes_[i]->A(A[i+2]*foot2);
	}
	exit.A(A[numNodes-1]*foot2);
	exit.Amax_ = A[numNodes-1]*foot2;


	hypro::systemModule::ModulePtr inlet(new hypro::FixedInlet("Inlet",&system.N1(),system.nodes_[0].get()));
	inlet->N1().Amax_ = M_PI*std::pow(0.263/2.0,2.0); //Maximum capture area (adapted conditions)
	inlet->N2().Name_ = "End of Intake";
	hypro::FixedInlet* inletC = (hypro::FixedInlet*)inlet.get();
	inletC->nodes_[0]->Name_ = "Intake Front";
	inletC->nodes_[0]->A(A[1]*foot2);
	hypro::ZitaInlet* inletC1 = (hypro::ZitaInlet*)inletC->modules_.front().get();
	//std::string directory(QUOTE(STXRescource));
	std::ifstream ifs (directory+"/InletProperties.dat");
	if(!ifs.good()){
		throw std::logic_error("Error: File InletProperties.dat not found in search path: "+directory+"/");
	}
	while(!ifs.eof() && ifs.good()){
		double zita_M, zita;
		ifs >> zita_M >> zita;
		inletC1->zita_[0].push_back(zita_M);
		inletC1->zita_[1].push_back(zita);
	}
	ifs.close();

	inletC->chokingFeedback_ = inlet;
	inletC->pdropmax_ = 200;
	inletC->unreduce();
	system.add(inlet);



	// First Ducting
	hypro::systemModule::ModulePtr ductone(new hypro::EffIsenDuct("First Ducting",system.nodes_[0].get(),system.nodes_[1].get()));
	ductone->N2().Name_ = "Pre-Injection";
	hypro::EffIsenDuct* ductoneC = (hypro::EffIsenDuct*)ductone.get();
	ductoneC->etap_ = 1;
	ductoneC->etaT_ = 1.0;
	ductoneC->chokingFeedback_ = inlet;
	system.add(ductone);


	// Injection
	hypro::systemModule::ModulePtr primary(new hypro::InjectionPhi<hypro::Friction>("primary",system.nodes_[1].get(),system.nodes_[2].get()));
	primary->N2().Name_ = "Injection Outlet";
	hypro::InjectionPhi<hypro::Friction>* primaryC = (hypro::InjectionPhi<hypro::Friction>*)primary.get();
	primaryC->PhiMax_ = 1.0;
	primaryC->Phi_ = 1.0;
	//List<double> Xf(5,0.0);
    std::vector<double> Xf = empty_thermos;
    Xf[53] = 1;
	primaryC->Rcoeff_[53] = 4;
	primaryC->Rcoeff_[1] = 71;
	primaryC->Xf_ = Xf;
	primaryC->Tf_ = 300;
	primaryC->mfrMax_ = fuelMfrMax;
	primaryC->chokingFeedback_ = inlet;
	primaryC->D_ = std::sqrt(primary->N2().A()*4/M_PI);
	primaryC->L_ = 1.0;
	primaryC->Cf_ = 4.0*primaryC->D_/(4.0*primaryC->L_);
	system.add(primary);


	// Secondary Ducting
	hypro::systemModule::ModulePtr ducttwo(new hypro::EffIsenDuct("Second Ducting",system.nodes_[2].get(),system.nodes_[3].get()));
	ducttwo->N2().Name_ = "Pre-Combustion";
	hypro::EffIsenDuct* ducttwoC = (hypro::EffIsenDuct*)ducttwo.get();
	ducttwoC->etap_ = 1;
	ducttwoC->etaT_ = 1.0;
	ducttwoC->chokingFeedback_ = inlet;
	system.add(ducttwo);


	// Combustor

	//systemModule::ModulePtr combustor(new NeutralLink("combustor",system.nodes_[3].get(),system.nodes_[4].get()));
	hypro::systemModule::ModulePtr combustor(new hypro::EffPhiComb<hypro::Friction>("combustor",system.nodes_[3].get(),system.nodes_[4].get()));
	combustor->N2().Name_ = "Combustor Outlet";
	combustor->chokingFeedback_ = primary;
	hypro::EffPhiComb<hypro::Friction>* combustorC = (hypro::EffPhiComb<hypro::Friction>*)combustor.get();
	std::ifstream ifs1 (directory+"/CombEff.dat");
	if(!ifs1.good()){
		throw std::logic_error("Error: File CombEff.dat not found in search path: "+directory+"/");
	}
	while(!ifs1.eof() && ifs1.good()){
		double eta_phi, eta;
		ifs1 >> eta_phi >> eta;
		combustorC->eta_[0].push_back(eta_phi);
		combustorC->eta_[1].push_back(eta);
	}
	ifs1.close();
	combustorC->Rcoeff_[53] = 4; // C12H23
	combustorC->Rcoeff_[3] = 71; // O2
	combustorC->Pcoeff_[5] = 46; // H2O
	combustorC->Pcoeff_[15] = 48; // CO2
	combustorC->iFuel_ = 53;
	combustorC->Cf_ = 0.0; //1.0;
	combustorC->D_ = 1.0;
	combustorC->L_ = 1.0;

	combustorC->chokingFeedback_ = inlet;
	system.add(combustor);



	// Nozzle Convergent

	hypro::systemModule::ModulePtr nozzleconv(new hypro::EffIsenDuct("NozzleConv",system.nodes_[4].get(),system.nodes_[5].get()));
	nozzleconv->N2().Name_ = "Nozzle Throat";
	hypro::EffIsenDuct* nozzleconvC = (hypro::EffIsenDuct*)nozzleconv.get();
	nozzleconvC->etap_ = 1.0; //0.96;
	nozzleconvC->etaT_ = 1.0;
	nozzleconvC->chokingFeedback_ = inlet;
	hypro::propModule::chokingMach_ = 1.0;
	system.add(nozzleconv);

	//Nozzle Divergent

	hypro::systemModule::ModulePtr nozzlediv(new hypro::EffIsenDuct("NozzleDiv",system.nodes_[5].get(),&system.N2()));
	hypro::EffIsenDuct* nozzledivC = (hypro::EffIsenDuct*)nozzlediv.get();
	nozzledivC->etap_ = 1.0; //0.96;
	nozzledivC->etaT_ = 1.0;
	nozzledivC->chokingFeedback_ = inlet;
	nozzledivC->choked_=true;
	system.add(nozzlediv);

	multiSystem.unreduce();

	hypro::systemModule& selMode(multiSystem.modes_[0]);
	selMode.verbosity_ = 1;
	selMode.calculate();
	Thrust = selMode.thrust();
	Mfr = selMode.propellantMfr();
	if(verbosity>0){
		system.printOut(std::cout);
	}

	std::vector<double> outTmp;
	outTmp.push_back(freeStream.A());
	outTmp.push_back(freeStream.Amax_);
	outTmp.push_back(freeStream.Amin_);
	outTmp.push_back(freeStream.getPress());
	outTmp.push_back(freeStream.getTemp());
	outTmp.push_back(freeStream.getU());
	for(std::size_t i=0; i<freeStream.X().size(); i++){
		outTmp.push_back(freeStream.X()[i]);
	}
	out.push_back(std::make_pair(freeStream.Name_,outTmp));

	hypro::Collection::OutType out1(selMode.out());
	for(std::size_t i=0; i<out1.size(); i++){
		out.push_back(out1[i]);
	}

	return new hypro::systemModule(system);
}
