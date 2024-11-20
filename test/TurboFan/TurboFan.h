/*!
// Created by Mark De Luca on 07/12/18.
//
// For any queries, contact m.de-luca.1@research.gla.ac.uk
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
 * along with HyPro.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <iostream>
#include <fstream>
#include <time.h>

#include "gasTurbine/Compressor.h"
#include "gasTurbine/Turbine.h"
#include "injectors/InjectionPhi.h"
#include "injectors/InjectionPhiNoTemp.h"

#include "combustion/EffComb.h"
#include "combustion/EffPhiComb.h"
#include "inletModels/FixedInlet.h"
#include "inletModels/Inlet.h"
#include "nozzleModels/AdaptedNozzle.h"
#include "nozzleModels/ConvNozzle.h"
#include "inletModels/AdaptedInlet.h"
#include "inletModels/AdaptedThroatInlet.h"
#include "core/Node.h"
#include "core/Collection.h"
#include "core/systemModule.h"
#include "solvers/EffIsenDuct.h"
#include "solvers/balanceMach.h"
#include "solvers/Friction.h"
#include "core/propModule.h"
#include "thermo/thermoState.h"
#include "NeutralLink.h"
#include "simpleSplitter.h"
#include "Mixer.h"

#include "cantera/thermo/ThermoPhase.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/Phase.h"
#include <memory>
#include "typedefs.h"
#include "offDesignData.h"


using namespace hypro;
std::pair<bool, outputData>
TurboFan (const double &p, const double &T, const double &Mach, const double &throttle, const int &verbosity, int &ondesign, double &CaptureArea, double &T04_, double &T045_, double &Aratio_HP,double &Aratio_LP, double &kH_HP, double &kH_LP) {

    outputData returned;

    //Node &freeStream =*new Node("/usr/share/cantera/data/gri30_highT.yaml", "gri30");
    Node &freeStream =*new Node("/home/andrew/GUHypro-tjet/test/TurboJet/gri30_highT_modified_jetA.yaml", "gri30");
    freeStream.Name_ = "Freestream";		// Sets the name of the freestream node
    
    Node &exit = *new Node(freeStream);
    exit.Name_ = "Nozzle exit to FS";		// Sets the name of the exit node

    // This creates a gas composition that approximates atmospheric conditions
    Cantera::compositionMap initialComposition;
    initialComposition.emplace("O2", 0.209);
    initialComposition.emplace("N2", 1 - 0.209);

    freeStream.setTPX(T, p, initialComposition);		// Sets the temperature, pressure, and composition of the freestream
    vector<double> empty_thermos(freeStream.X().size());	// Creates an empty vector, and sets the size to the no. species within the gas mixture
    freeStream.setU(Mach * freeStream.geta());		// Calculates the flow velocity of the freestream

    systemModule tfan("system", &freeStream, &exit);		// Creates an object which models a system that involves the inlet and exit nozzle flows
    tfan.verbosity_ = verbosity;				// Sets the verbosity level of the tjet systemModule

    const int numNodes = 7;			// Sets the no. nodes within the system
    static vector<double> A(numNodes) ;	// Creates an empty vector set to contain the no. elements dictated by variable numNodes

	const double AA = 5.725;	// THIS VALUE SETS THE AREA RATIO
//	const double AA = 1;

//	const double AAprime = 5.4394;
    A = {AA, AA, AA, AA, AA, AA, 0.52};

    const int numInternalNodes = numNodes - 2;

    // Sets the area of the freestream node	
    if (CaptureArea){				// If non-zero
        freeStream.A(A[0]*CaptureArea);
    }
    else{
        freeStream.A(A[0]);
    }

    tfan.initNodes(numInternalNodes, freeStream);
    for (int i = 0; i < numInternalNodes; i++) {		// Sets the area of the internal nodes
    	tfan.nodes_[i]->A(A[i + 1]);
    }
    
    exit.A(A[numNodes - 1]);		// Sets the area of the exit node
//    exit.Amax_=AA;


    const int numShafts = 1;
    tfan.initShaft(numShafts);
    
    double ERR = 100;
    double pi_fan = 1.5;
    double pif_max = 2.5;
    double pif_min = 1.01;


/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Set up of Engine Configuration
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/

    //double Tfan
    //systemModule::ModulePtr fan(new Fan("fan",))

    double Tfan ();

    systemModule::ModulePtr inlet(new Inlet("inlet", &tfan.N1(), tfan.nodes_[0].get()));
    Inlet *inletC = (Inlet *) inlet.get();
    inlet->N2().Name_ = "LP Compressor Inlet";
    tfan.add(inlet);

//    systemModule::ModulePtr inlet(new AdaptedInlet("inlet", &tfan.N1(), tfan.nodes_[0].get()));	// THIS LIMITS M_inf TO 0.9
//    AdaptedInlet *inletC = (AdaptedInlet *) inlet.get();
//    inlet->N2().Name_ = "LP Compressor Inlet";
//    tfan.add(inlet);


    systemModule::ModulePtr lpcompressor(new Compressor("lpcompressor", tfan.nodes_[0].get(), tfan.nodes_[1].get(), tfan.mech_[0].get()));
    Compressor *lpcompressorC = (Compressor *) lpcompressor.get();
    lpcompressorC->effc_ = 0.8782;
    lpcompressorC->pi_c_ = 4.1;
    lpcompressorC->checkdesignpoint_ = ondesign;
    tfan.mech_[0]->eff_mech_ = 0.99;		// NOT IN TURBOJET.H
    if (ondesign == 0){
    	lpcompressorC->setT04(T045_);
    	lpcompressorC->kH_ = kH_LP;
        tfan.mech_[0]->Cp2_ =  tfan.nodes_[0]->Cp();
    	tfan.mech_[0]->Cp_ =  tfan.nodes_[0]->Cp();
    	tfan.mech_[0]->mfr_ =  tfan.nodes_[0]->mfr();
    	tfan.mech_[0]->kH_ = kH_LP;
    	tfan.mech_[0]->Ht_ = 184286*freeStream.mfr();
		}
//    lpcompressorC->chokingFeedback_ = inlet;
    tfan.add(lpcompressor);						// Adds the module 'compressor' to the tjet system

    systemModule::ModulePtr hpcompressor(new Compressor("hpcompressor", tfan.nodes_[1].get(), tfan.nodes_[2].get(), tfan.mech_[1].get()));
    Compressor *hpcompressorC = (Compressor *) hpcompressor.get();
    hpcompressorC->effc_ = 0.8402;
    hpcompressorC->pi_c_ = 2.9;
    hpcompressorC->checkdesignpoint_ = ondesign;
    tfan.mech_[1]->eff_mech_ = 0.99;
    if (ondesign == 0){
    	hpcompressorC->setT04(T04_);
    	hpcompressorC->kH_ = kH_HP;
    	tfan.mech_[1]->Cp2_ =  tfan.nodes_[0]->Cp();
    	tfan.mech_[1]->Cp_ =  tfan.nodes_[0]->Cp();
    	tfan.mech_[1]->mfr_ =  tfan.nodes_[0]->mfr();
    	tfan.mech_[1]->kH_ = kH_HP;
    	tfan.mech_[1]->Ht_ = 267345*freeStream.mfr();
		}
    tfan.add(hpcompressor);
    
    
    systemModule::ModulePtr injection(new InjectionPhiNoTemp<balanceMach>("injection", tfan.nodes_[2].get(),  tfan.nodes_[3].get()));
    InjectionPhiNoTemp<balanceMach> *injectionC = (InjectionPhiNoTemp<balanceMach> *) injection.get();
    injectionC->Phi_ = 1.0 * throttle;
    injectionC->PhiMax_ = 1.0;
    injectionC->chokingFeedback_ = injection;
    injectionC->chokingMach_ = 0.98;
    injectionC->Tf_ = 854;
    vector<double> Xf = empty_thermos;
//    Xf[53] = 1;
//	injectionC->Rcoeff_[3] = 71;
//	injectionC->Rcoeff_[53] = 4;
	Xf[5] = 1;
	injectionC->Rcoeff_[0] = 71;
	injectionC->Rcoeff_[5] = 4;
    injectionC->Xf_ = Xf;
    tfan.add(injection);


    systemModule::ModulePtr combustion(new EffComb<balanceMach>("combustor", tfan.nodes_[3].get(), tfan.nodes_[4].get()));
    EffComb<balanceMach> *combustionC = (EffComb<balanceMach> *) combustion.get();
//	combustionC->Rcoeff_[53] = 4;
//	combustionC->Rcoeff_[3]  = 71;
//	combustionC->Pcoeff_[15] = 48;
//	combustionC->Pcoeff_[5]  = 46;
	combustionC->Rcoeff_[5] = 4;
	combustionC->Rcoeff_[0]  = 71;
	combustionC->Pcoeff_[2] = 48;
	combustionC->Pcoeff_[1]  = 46;
	combustionC->eta_ = 1.0;
	combustionC->chokingFeedback_ = injection;
	tfan.add(combustion);			// Adds the module 'combustion' to the tjet system


    systemModule::ModulePtr hpturbine(new Turbine("hpturbine", tfan.nodes_[4].get(), tfan.nodes_[5].get(), tfan.mech_[1].get()));
    Turbine *hpturbineC = (Turbine *) hpturbine.get();
    hpturbineC->efft_ = 0.8785;
//    hpturbineC->efft_ = 0.91325;
    hpturbineC->checkdesignpoint_ = ondesign;
	if (ondesign == 0){
    		hpturbineC->Aratio_ = Aratio_HP;
    		hpturbineC->kH_ = kH_HP;
    	    hpturbineC->efft_ = 0.8785;
	}
    tfan.add(hpturbine);


   systemModule::ModulePtr lpturbine(new Turbine("lpturbine", tfan.nodes_[5].get(),  tfan.nodes_[6].get(), tfan.mech_[0].get()));
    Turbine *lpturbineC = (Turbine *) lpturbine.get();
//    lpturbineC->efft_ = 0.922; //Nozzle throat area is 0.7725
    lpturbineC->efft_ = 0.8890; //Nozzle throat area is 0.79981
    lpturbineC->checkdesignpoint_ = ondesign;
	if (ondesign == 0){
    		lpturbineC->Aratio_ = Aratio_LP;
    		lpturbineC->kH_ = kH_LP;
    	    lpturbineC->efft_ = 0.8890;
	}
    tfan.add(lpturbine);
    
    
    systemModule::ModulePtr nozzle(new AdaptedNozzle<EffIsenDuct>("nozzle", tfan.nodes_[6].get(), &tfan.N2()));
    AdaptedNozzle<EffIsenDuct> *nozzleC = (AdaptedNozzle<EffIsenDuct> *) nozzle.get();	// NozzleC is a pointer to the AdaptedNozzle template class which is stored inside the pointer to the 'nozzle' object - allows directed access of AdaptedNozzle template class.
  	nozzleC->etap_ = 1.0;
  	nozzleC->chokingMach_ = 1.0;			// Could this be changed - only works for up to M = 1.0
  	nozzleC->freeStream_ =  &inlet->N1();		// Sets the freestream conditions for the nozzle
  	tfan.add(nozzle);


//    systemModule::ModulePtr nozzle(new AdaptedNozzle<EffIsenDuct>("nozzle", tfan.nodes_[14].get(), &tfan.N2()));
//    AdaptedNozzle<EffIsenDuct> *nozzleC = (AdaptedNozzle<EffIsenDuct> *) nozzle.get();
//  	nozzleC->etap_ = 1.0;
//  	nozzleC->chokingMach_ = 0.99;
//  	nozzleC->freeStream_ =  &freeStream;
//  	tfan.add(nozzle);


  	tfan.calculate();
  	returned.sys = new systemModule(tfan);
    returned.thrust =tfan.thrust();
    returned.Mfr =tfan.propellantMfr();

  	if (!((CaptureArea == false) && (ondesign == 0))) {

  		tfan.printOut(std::cout);
  	  	std::cout << "\n" << "Thrust is: " << tfan.thrust() << " N" << "\n";

	}


    returned.T04 = combustion->N2().getT0();
   	returned.Aratio_HP = hpturbineC->Aratio_;
   	returned.kH_HP = hpturbineC->kH_;
    returned.T045 = hpturbine->N2().getT0();
    returned.Aratio_LP = lpturbineC->Aratio_;
    returned.kH_LP = lpturbineC->kH_;
}





















