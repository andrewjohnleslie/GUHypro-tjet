/*!
 *  \brief      Validation test case based on the Hyperion propulsion system
 *  \details    Hyperion is a space plane vehicle developed at the GeorgiaTech.
 *              Details of the validation are available in the
 *              <a href="http://oleg.lib.strath.ac.uk/R/?func=dbin-jump-full&object_id=26892">PhD thesis of Alessandro Mogavero</a>
 *
 *              @param[in] p free stream static pressure
 *              @param[in] T free stream static temperature
 *              @param[in] Mach free stream Mach number
 *              @param[in] Throttle
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

//#include "core/Node.h"
#include "Node.h"

#include "core/systemModule.h"
#include "core/MultiModeMachP0.h"
#include "solvers/isentropicDuct.h"
#include "solvers/EffIsenDuct.h"
#include "solvers/Friction.h"
//#include "ISAModel.h"
#include "thermo/thermoState.h"
#include "NeutralLink.h"
 
#include "Wedge.h"
#include "inletModels/AdaptedThroatInlet.h"
#include "inletModels/ScramjetInlet.h"
#include "inletModels/SupersonicInlet.h"
#include <AdaptedInlet.h>
#include "Rocket.h"
#include "Mixer.h"
#include "injectors/InjectionPhi.h"
#include "combustion/EffComb.h"
#include "combustion/EquilCombustion.h"
#include "nozzleModels/AdaptedNozzle.h"

#include "cantera/thermo/ThermoPhase.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/Phase.h"
#include <memory>
//#include "fmt/format.h"
#include "returnData.h"
using namespace hypro;

//	systemModule* Hyperion(const double& p, const double& T, const double& Mach, const double& throttle, const int mode, const double& verbosity,
//						   double& Thrust, vector<double>& Mfr, Collection::OutType& out, double& dynamic_pressure, vector<double> &emissions, double &total_flow){
std::pair<bool, outputData>
Chiara(const double &p, const double &T, const double &Mach, const double &throttle, const int mode,
         const int &verbosity) {
/** Modes:
		-1 - Automatic switching
		0 - Ejector Mode
		1 - Ramjet Mode
		2 - Scramjet Mode
		3 - Rocket Mode
**/
    outputData returned;
    double foot = 0.3048; //conversion feet to meters
    double lbm = 0.45359; //conversion libre to kg

    Node freeStream("/usr/share/cantera/data/gri30_highT.yaml", "gri30");
    double pStream = p;
    if (p < 1) {
        pStream = 1;
    }
    //freeStream.p_ = pStream;
    //freeStream.T_ = T;
///////////////////////////////////////////////////////////
    Cantera::compositionMap initialComposition;
    initialComposition.emplace("O2", 0.209);
    initialComposition.emplace("N2", 1 - 0.209);

    freeStream.setTPX(T, pStream, initialComposition);
    vector<double> empty_thermos(freeStream.X().size());

//    // For ramjet, if mach < 1 a nan is the returned thrust. between mach 1 and
//    // ~1.5 the returned thrust is negative. If in ramjet mode change to mach 1
//    // if mach < 1
//    double MachStream = Mach;
//    if (mode == 1) {
//        if (Mach < 1) {
//            MachStream = 1;
//        }
//    }

    freeStream.setU(Mach * freeStream.geta());
    freeStream.Name_ = "Free Stream";

    thermoState::warning_ = false;

//The two following lines mean that the default name of all the nodes is set to "Pre Intake"
    Node &freeStreamMod = *new Node("/usr/share/cantera/data/gri30_highT.yaml", "gri30");     
    freeStreamMod.Name_ = "Pre Intake";
    Node &exit = *new Node("/usr/share/cantera/data/gri30_highT.yaml", "gri30");
    exit.Name_ = "Nozzle Exit";

////////////////////////////////////////////////

    Wedge forebody("forebody", &freeStream, &freeStreamMod);
    forebody.delta_ = 8.0;
    forebody.calculate();

    MultiModeMachP0 multiSystem("system", &freeStreamMod, &exit);
    multiSystem.verbosity_ = verbosity;
    systemModule &system(multiSystem.add("Ejector Mode"));
    system.verbosity_ = verbosity;

//////////////////////////////////////////////////////

 //Read cross sectional areas from file
    const int numNodes = 12;
    double A[numNodes];

    ifstream inFile;
    inFile.open("/home/mark/git/GUHypro/build/test/Chiara/Input_A.txt");
  if(!inFile){
	 cout << "Error in opening file" << endl;
  }else{
	  for (int i=0; i<numNodes;i++){
		  inFile >> A[i];
		//cout << A[i] << endl;
	  }
	 inFile.close();
}

//Nodes                 N1_    I0    I1       0     1       2     3      4        5      6     7     N2_95
//
//double A[numNodes] =   {27.0, 27.0, 0.25*27, 8.24, 11.25, 11.25, 0.0, 11.25-8.24, 22.5, 22.5, 22.5, 44.3};
//
//    /*name = {'Pre Intake','Intake','Throat','Pinch Point','Primary','Mixer End',...
//    'Dummy','Rocket Outlet','End Diffuser','Injection','End Chamber','Nozzle'};*/
    double foot2 = std::pow(foot, 2.0);
    //std::cout << A[0]*foot2 << std::endl;
    // std::cout << A[1]*foot2 << std::endl;
    freeStreamMod.A(A[0] * foot2);
    system.initNodes(8, freeStreamMod);
    for (int i = 0; i < 8; i++) {
        system.nodes_[i]->A(A[i + 3] * foot2);
    }
    exit.A(A[numNodes - 1] * foot2);
    exit.Amax_ = A[numNodes - 1] * foot2;

    ///////////////////////////////////////////

    //======================== Ejector Mode =============================
    //Inlet
    systemModule::ModulePtr inlet(new AdaptedThroatInlet("Inlet", &system.N1(), system.nodes_[0].get()));
    inlet->N2().Name_ = "Pinch Point";
    AdaptedThroatInlet *inletC = (AdaptedThroatInlet *) inlet.get();
    inletC->nodes_[0]->A(A[1] * foot2 * throttle);
    inletC->nodes_[0]->Amax_ = inletC->nodes_[0]->A();
    inletC->nodes_[0]->Name_ = "Intake";
    inletC->nodes_[1]->A(A[2] * foot2 * throttle);
    inletC->nodes_[1]->Amax_ = inletC->nodes_[1]->A();
    inletC->nodes_[1]->Name_ = "Throat";
    inlet->verbosity_ = 0;
    inlet->chokingFeedback_ = inlet;
    system.add(inlet);



    // Primary - ejector
    systemModule::ModulePtr postPinch(new isentropicDuct(
            "Post pinch point", &inlet->N2(), system.nodes_[1].get()));
    system.add(postPinch);

//the node system.nodes_[3] is used as the "upstream" start of the rocket.Its actual values have no effecto on the rocket, its just that the rocket "module" needs a start point by definition.

    systemModule::ModulePtr rocket(new Rocket("Rocket", system.nodes_[3].get(), system.nodes_[4].get()));
    rocket->N1().Name_ ="Rocket Inlet";
    rocket->N2().Name_ = "Rocket exit";
    Rocket *rocketC = (Rocket *) rocket.get();
    vector<double> Xr = empty_thermos;
    Xr[5] = 1;
    rocketC->N2().X(Xr);
    rocketC->Me_ = 6.0;
    rocketC->mfrMax_ = 216 * lbm * throttle; //TODO Add throttle in here
    rocketC->Isp_ = 462;
    rocketC->unreduce();
    system.add(rocket);

    systemModule::ModulePtr primary(new Mixer("primary",
                                              &postPinch->N2(), system.nodes_[2].get(), &rocket->N2()));
    primary->N1().Name_ = "Primary In";
    primary->N2().Name_ = "Primary Out";
    Mixer *primaryC = (Mixer *) primary.get();
    primaryC->verbosity_ = 0;
    primaryC->chokingFeedback_ = inlet;
    primaryC->eff_ = 0.9;
    system.add(primary);

    // Secondary - Post combustor
    systemModule::ModulePtr diffuser(new EffIsenDuct("diffuser", system.nodes_[2].get(), system.nodes_[5].get()));
    diffuser->N2().Name_ = "DiffuserEnd";   //levo lo spazio tra Chamber ed End per poter lavorare meglio sulle celle quando importo i dati in Matlab
    EffIsenDuct *diffuserC = (EffIsenDuct *) diffuser.get();
    diffuserC->etap_ = 0.92;
    diffuserC->etaT_ = 1.0;
    system.add(diffuser);

    systemModule::ModulePtr secondary(
            new InjectionPhi<balanceMach>("secondary", system.nodes_[5].get(), system.nodes_[6].get()));
    secondary->N2().Name_ = "Injection";
    InjectionPhi<balanceMach> *secondaryC = (InjectionPhi<balanceMach> *) secondary.get();
    secondaryC->PhiMax_ = 1 * throttle; //TODO Add throttle for injector here
    secondaryC->Phi_ = 1 * throttle;
    //secondaryC->PhiMax_ = 1; //TODO Add throttle for injector here
    //secondaryC->Phi_ = 1;
    //secondaryC->mfrMax_ = *throttle;
    vector<double> Xf = empty_thermos;
    Xf[0] = 1;
    secondaryC->Xf_ = Xf;
    secondaryC->Tf_ = 20;
    secondaryC->Rcoeff_[0] = 2;
    secondaryC->Rcoeff_[3] = 1;
    secondaryC->verbosity_ = verbosity;
    secondary->chokingFeedback_ = secondary;
    secondaryC->chokingMach_ = 0.9;
    system.add(secondary);

	systemModule::ModulePtr combustor(new EffComb<Friction>("combustor",system.nodes_[6].get(),system.nodes_[7].get()));
	combustor->N2().Name_ = "ChamberEnd";  //stesso discorso valido per Chamber e End
	EffComb<Friction>* combustorC = (EffComb<Friction>*)combustor.get();
	combustorC->Rcoeff_[0] = 2;
	combustorC->Rcoeff_[3] = 1;
	combustorC->Pcoeff_[5] = 2;
	combustorC->eta_ = 0.8;
	combustorC->Cf_ = 0.01;
	combustorC->D_ = 1.631;
	combustorC->L_ = 1.0;
	combustorC->chokingFeedback_ = secondary;
	system.add(combustor);


//     systemModule::ModulePtr combustor(
//                 new EquilCombustion<Friction>("combustor", system.nodes_[6].get(), system.nodes_[7].get()));
//     combustor->N2().Name_ = "Combustor Outlet";
//     EquilCombustion<Friction> *combustorC = (EquilCombustion<Friction> *) combustor.get();
//     combustorC->verbosity_ = verbosity - 1;
//     combustorC->chokingFeedback_ = secondary;
//     combustorC->Cf_ = 0.01;
//	 combustorC->Cf_ = 0.00;
//     combustorC->D_ = 1.631;
//     combustorC->L_ = 1.0;
//     system.add(combustor);

//	systemModule::ModulePtr combustorw

//    systemModule::ModulePtr combustor(
//            new EquilCombustion<balanceMach>("combustor", system.nodes_[6].get(), system.nodes_[7].get()));
//    combustor->N2().Name_ = "Combustor Outlet";
//    EquilCombustion<balanceMach> *combustorC = (EquilCombustion<balanceMach> *) combustor.get();
//    combustorC->verbosity_ = verbosity - 1;
//    combustorC->chokingFeedback_ = secondary;
////        combustorC->Cf_ = 0.01;
////        combustorC->Cf_ = 0.00;
////        combustorC->D_ = 1.631;
////        combustorC->L_ = 1.0;
//    system.add(combustor);

    // Nozzle
    systemModule::ModulePtr nozzle(new AdaptedNozzle<EffIsenDuct>("nozzle", system.nodes_[7].get(), &system.N2()));
    AdaptedNozzle<EffIsenDuct> *nozzleC = (AdaptedNozzle<EffIsenDuct> *) nozzle.get();
    nozzleC->choked_ = true;
    /* in case of inlet spilling (especially in supersonic) better to adapt the
        nozzle to the inlet conditions.*/
    nozzleC->freeStream_ = &inlet->N1();
    nozzleC->etap_ = 1.0;
    system.add(nozzle);


    // ========================= Ramjet mode =========================
    systemModule &ramjet(multiSystem.add("Ramjet Mode", 0));
    const systemModule::ModulePtr inlet1(ramjet.exchange<AdaptedThroatInlet>(inlet, "Inlet"));
    AdaptedThroatInlet *inlet1C = (AdaptedThroatInlet *) inlet1.get();
    inlet1C->nodes_[0]->A(A[1] * foot2);
    inlet1C->nodes_[0]->Amax_ = inlet1C->nodes_[0]->A();
    inlet1C->nodes_[0]->Name_ = "Intake";
    inlet1C->nodes_[1]->A(1.055 * A[2] * foot2);
    inlet1C->nodes_[1]->Amax_ = inlet1C->nodes_[1]->A();
    inlet1C->nodes_[1]->Name_ = "Throat";
    inlet->verbosity_ = 0;
    inlet->chokingFeedback_ = inlet;

    const systemModule::ModulePtr primary1(ramjet.exchange<EffIsenDuct>(primary, "primary ramjet"));
    primary1->chokingFeedback_ = inlet1;
    EffIsenDuct *primary1C = (EffIsenDuct *) primary1.get();
    primary1C->etap_ = 0.725;

    ramjet.remove(rocket);


    //========================= Scramjet mode =========================
    systemModule &scramjet(multiSystem.add("Scramjet Mode", 0));
    scramjet.verbosity_ = 1;
    /*const systemModule::ModulePtr inlet2(scramjet.exchange<EffIsenDuct>(inlet,"Inlet scramjet"));
    EffIsenDuct* inlet2C = (EffIsenDuct*)inlet2.get();
    inlet2C->etap_ = 0.5;
    inlet2C->verbosity_ =1;
    inlet2C->choked_ = false;*/

    const systemModule::ModulePtr inlet2(scramjet.exchange<ScramjetInlet>(inlet, "Inlet"));
    ScramjetInlet *inlet2C = (ScramjetInlet *) inlet2.get();
    inlet2C->verbosity_ = 1;

    const systemModule::ModulePtr primary3(scramjet.exchange<EffIsenDuct>(primary, "primary scramjet"));
    EffIsenDuct *primary3C = (EffIsenDuct *) primary3.get();
    primary3C->etap_ = 1;

    /*const systemModule::ModulePtr secondary3(scramjet.exchange<InjectionPhi<Friction>>(secondary,"secondary scramjet"));
    InjectionPhi<Friction>* secondary3C = (InjectionPhi<Friction>*)secondary3.get();
    secondary3C->PhiMax_ = 1*throttle;
    secondary3C->Phi_ = 1*throttle;
    //secondaryC->mfrMax_ = *throttle;
    //List<double> Xf(4,0.0);
    //Xf[0] = 1;
//	secondary3C->Xf_ = Xf;
//	secondary3C->Tf_ = 20;
//	secondary3C->Rcoeff_[0] = 2;
//	secondary3C->Rcoeff_[3] = 1;
    secondary3C->Cf_ = 0.035;
    secondary3C->D_ = 4*0.026*1/(2*1);
    secondary3C->L_ = 0.25;
    secondary3C->chokingMach_ = 1.1;

//    vector<double> Xf = empty_thermos;
//    Xf[0] = 1;
    secondary3C->Xf_ = Xf;
    secondary3C->Tf_ = 20;
    secondary3C->Rcoeff_[0] = 2;
    secondary3C->Rcoeff_[3] = 1;
    secondary3->chokingFeedback_ = secondary3;*/

    /*const systemModule::ModulePtr combustor2(scramjet.exchange<EffComb<Friction>>(combustor,"combustor scramjet"));
    EffComb<Friction>* combustor2C = (EffComb<Friction>*)combustor2.get();
    combustor2C->N2().Name_ = "Chamber End";
    combustor2C->chokingFeedback_ = secondary;
    combustor2C->Rcoeff_[0] = 2;
    combustor2C->Rcoeff_[3] = 1;
    combustor2C->Pcoeff_[5] = 2;
    combustor2C->eta_ = 0.8;
    combustor2C->Cf_ = 0.035;
    combustor2C->D_ = 1.631;
    combustor2C->L_ = 1.0;

    const systemModule::ModulePtr combustor2(scramjet.exchange<EquilCombustion<Friction>>(combustor,"combustor scramjet"));
    EquilCombustion<Friction>* combustor2C = (EquilCombustion<Friction>*)combustor2.get();
    combustor2C->N2().Name_ = "Combustor Outlet";
    combustor2C->verbosity_ = verbosity - 1;
    combustor2C->chokingFeedback_ = secondary;
    combustor2C->Cf_ = 0.01;
    combustor2C->D_ = 1.631;
    combustor2C->L_ = 1.0;
    system.add(combustor);*/

    scramjet.remove(rocket);

    //========================= Pure Rocket mode
    systemModule &rocketMode(multiSystem.add("Rocket Mode", 0));
    const systemModule::ModulePtr primary2(rocketMode.exchange<isentropicDuct>(primary, "primary"));
    primary2->N1(rocket->N2());
    const systemModule::ModulePtr secondary2(rocketMode.exchange<NeutralLink>(secondary, "secondary"));
    secondary2->N2(combustor->N2());
    rocketMode.remove(inlet);
    rocketMode.remove(postPinch);
    rocketMode.remove(combustor);
    Node rocketDummyN1(rocketMode.N1());
    rocketMode.N1(rocketDummyN1);
    rocketDummyN1.setU(0.0);



    multiSystem.rangesM_[0][0] = 0.0;
    multiSystem.rangesM_[0][1] = 3.0;
    multiSystem.rangesM_[1][0] = 3.0;
    multiSystem.rangesM_[1][1] = 7.0;
    multiSystem.rangesP0_[0][0] = 1e4;
    multiSystem.rangesP0_[1][0] = 5e5;

    if (mode == -1) {
        // Automatic mode switching
        multiSystem.calculate();
        returned.thrust = multiSystem.thrust();
        returned.Mfr = multiSystem.propellantMfr();
//		Thrust = multiSystem.thrust();
//		Mfr = multiSystem.propellantMfr();
    } else if (mode >= 0) {
        systemModule &selMode(multiSystem.modes_[mode]);
        selMode.verbosity_ = verbosity;
        selMode.calculate();

        returned.thrust = selMode.thrust();
        returned.Mfr = selMode.propellantMfr();
        returned.dynamic_pressure = (0.5 * freeStream.rho() * std::pow(freeStream.getU(), 2));

        returned.total_flow = system.N2().rho() * system.N2().getU() * system.N2().A();
        returned.emissions.push_back(system.N2().X("H2O"));
        returned.emissions.push_back(system.N2().X("N2"));
        returned.emissions.push_back(system.N2().X("O2"));
        returned.emissions.push_back(system.N2().X("OH"));
        returned.emissions.push_back(system.N2().X("NO"));
        returned.emissions.push_back(system.N2().X("NO2"));
        Collection::OutType out1(selMode.out());


        for (std::size_t i = 0; i < out1.size(); i++) {
            returned.out.push_back(out1[i]);
        }

		selMode.printOut(cout);
		selMode.saveJson("/home/mark/git/GUHypro/build/test/Chiara/Chiara.txt");
    }
    returned.sys = new systemModule(system);
    return {true, returned};
}
