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


#include "gasTurbine/Compressor.h"	//DONE
#include "gasTurbine/Turbine.h"	//DONE
#include "injectors/InjectionPhi.h"
#include "combustion/EffComb.h"
#include "combustion/EffPhiComb.h"
#include "inletModels/FixedInlet.h"
#include "inletModels/AdaptedInlet.h"
#include "nozzleModels/AdaptedNozzle.h"
#include "nozzleModels/ConvNozzle.h"


#include "core/Node.h"
#include "core/Collection.h"
#include "core/systemModule.h"
#include "solvers/EffIsenDuct.h"
#include "solvers/balanceMach.h"
#include "solvers/Friction.h"
#include "core/propModule.h"
#include "thermo/thermoState.h"
#include "NeutralLink.h"


#include "cantera/thermo/ThermoPhase.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/Phase.h"
#include <memory>
#include "typedefs.h"
#include "offDesignData.h"




using namespace hypro;
std::pair<bool, outputData>
TurboJet(const double &p, const double &T, const double &Mach, const double &throttle, const int &verbosity, int &ondesign, double &CaptureArea, double &T04_, double &Aratio_, double &kH_) {

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

    systemModule tjet("system", &freeStream, &exit);		// Creates an object which models a system that involves the inlet and exit nozzle flows
    tjet.verbosity_ = verbosity;				// Sets the verbosity level of the tjet systemModule

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

    tjet.initNodes(numInternalNodes, freeStream);
    for (int i = 0; i < numInternalNodes; i++) {		// Sets the area of the internal nodes
    	tjet.nodes_[i]->A(A[i + 1]);
    }
    
    exit.A(A[numNodes - 1]);		// Sets the area of the exit node
//    exit.Amax_=AA;


    const int numShafts = 1;
    tjet.initShaft(numShafts);


/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Set up of Engine Configuration
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/

    systemModule::ModulePtr inlet(new EffIsenDuct("inlet", &tjet.N1(), tjet.nodes_[0].get()));
    inlet->N2().Name_ = "Compressor Inlet";			// Accesses the second node of the 'inlet' and sets its name.
    EffIsenDuct *inletC = (EffIsenDuct *) inlet.get();	// InletC is a pointer to the EffIsenDuct class which is stored inside the pointer to the 'inlet' object 
    									// Allows direct access to members and methods of the EffIsenDuct class.
    inletC->etap_ = 1.0;
	inletC->etaT_ = 1.0;
	inletC->choked_ = false;
	tjet.add(inlet);					// Adds the module 'inlet' to the tjet system.

//    systemModule::ModulePtr inlet(new Inlet("inlet", &tjet.N1(), tjet.nodes_[0].get()));
//    Inlet *inletC = (Inlet *) inlet.get();
//    inlet->N2().Name_ = "LP Compressor Inlet";
//    tjet.add(inlet);

//    systemModule::ModulePtr inlet(new AdaptedInlet("inlet", &tjet.N1(), tjet.nodes_[0].get()));	// THIS LIMITS M_inf TO 0.9
//    AdaptedInlet *inletC = (AdaptedInlet *) inlet.get();
//    inlet->N2().Name_ = "LP Compressor Inlet";
//    tjet.add(inlet);


    systemModule::ModulePtr compressor(new Compressor("compressor", tjet.nodes_[0].get(), tjet.nodes_[1].get(), tjet.mech_[0].get()));
    compressor->N2().Name_ = "Compressor Exit";			// Accesses the second node of the 'compressor' and sets its name.
    Compressor *compressorC = (Compressor *) compressor.get();	// CompressorC is a pointer to the Compressor class which is stored inside the pointer to the 'compressor' object - allows directed access of Compressor class.
    compressorC->effc_ = 0.9;
    compressorC->pi_c_ = 5.5;
//    compressorC->pi_c_ = 5.0323585;
    compressorC->checkdesignpoint_ = ondesign;
    if (!ondesign){
			compressorC->setT04(T04_);
			tjet.mech_[0]->Cp2_ =  tjet.nodes_[0]->Cp();	// Brings the variables within the scope of the mechanicalLink class
			tjet.mech_[0]->Cp_ =  tjet.nodes_[0]->Cp();
			tjet.mech_[0]->mfr_ =  tjet.nodes_[0]->mfr();
			tjet.mech_[0]->kH_ = kH_;
			tjet.mech_[0]->Ht_ = 258307*freeStream.mfr();
		}
	tjet.add(compressor);						// Adds the module 'compressor' to the tjet system


    systemModule::ModulePtr injection(new InjectionPhi<balanceMach>("injection", tjet.nodes_[1].get(),  tjet.nodes_[2].get()));
    injection->N2().Name_ = "Combustion Inlet";						// Accesses the second node of the 'combustion' and sets its name.
    InjectionPhi<balanceMach> *injectionC = (InjectionPhi<balanceMach> *) injection.get();	// InjectionC is a pointer to the InjectionPhi template class which is stored inside the pointer to the 'injection' object - allows directed access of InjectionPhi template class.
    injectionC->Phi_ = 1.0 * throttle;
    injectionC->PhiMax_ = 1.0;
    injectionC->chokingFeedback_ = injection;
    injectionC->chokingMach_ = 0.98;
//    if (ondesign == 0){
//    	injectionC->Tf_ = 492.08;
//    }
//    else{
//    	injectionC->Tf_ = 479.18;
//    }
    if (!ondesign){
    	injectionC->Tf_ = 492.08;			// Sets the fuel injection temperature
    }
    else{
    	injectionC->Tf_ = 492.08;
    }
    vector<double> Xf = empty_thermos;
    Xf[53] = 1;
	injectionC->Rcoeff_[3] = 71;			// Sets the reaction coefficient for specifies 
	injectionC->Rcoeff_[53] = 4;
	injectionC->Xf_ = Xf;
    tjet.add(injection);				// Adds the module 'injection' to the tjet system



    systemModule::ModulePtr combustion(new EffComb<balanceMach>("combustor", tjet.nodes_[2].get(),  tjet.nodes_[3].get()));
    combustion->N2().Name_ = "Turbine Inlet";						// Accesses the second node of the 'combustion' and sets its name.
    EffComb<balanceMach> *combustionC = (EffComb<balanceMach> *) combustion.get();	// CombustionC is a pointer to the EffComb template class which is stored inside the pointer to the 'combustion' object - allows directed access of EffComb template class.
	combustionC->Rcoeff_[53] = 4;			// Sets the reaction coefficient for specifies 
	combustionC->Rcoeff_[3]  = 71;
	combustionC->Pcoeff_[15] = 48;
	combustionC->Pcoeff_[5]  = 46;
	combustionC->eta_ = 1.0;
	combustionC->chokingFeedback_ = injection;
	tjet.add(combustion);				// Adds the module 'combustion' to the tjet system




    /*systemModule::ModulePtr injection(new InjectionPhi<Friction>("injection", tjet.nodes_[1].get(),  tjet.nodes_[2].get()));
    injection->N2().Name_ = "Combustion Inlet";
    InjectionPhi<Friction> *injectionC = (InjectionPhi<Friction> *) injection.get();
    injectionC->Phi_ = 1.0 * throttle;
    injectionC->PhiMax_ = 1.0;
    injectionC->chokingFeedback_ = injection;
    injectionC->chokingMach_ = 0.98;
    injectionC->Tf_ = 300;
    vector<double> Xf = empty_thermos;
    Xf[0] = 1;
  	injectionC->Rcoeff_[0] = 2;
  	injectionC->Rcoeff_[3] = 1;
  	injectionC->Cf_ = 0.0035; //1.0;
  	injectionC->D_  = 4*(2*AA)*1/(2*1);
	injectionC->L_  = 0.25;
    injectionC->Xf_ = Xf;
    tjet.add(injection);



    systemModule::ModulePtr combustion(new EffComb<Friction>("combustor", tjet.nodes_[2].get(),  tjet.nodes_[3].get()));
    combustion->N2().Name_ = "Turbine Inlet";
    EffComb<Friction> *combustionC = (EffComb<Friction> *) combustion.get();
  	combustionC->eta_ = 1.0;
  	combustionC->Rcoeff_[0] = 2;
  	combustionC->Rcoeff_[3] = 1;
  	combustionC->Pcoeff_[5] = 2;
	combustionC->Cf_ = 0.0035; //1.0;
	combustionC->D_  = 4*(2*AA)*1/(2*1);
	combustionC->L_  = 0.25;
	tjet.add(combustion);


    //Comment in if wishing to use Hydrogen as fuel type instead of Jet-A (C12H23) */


    systemModule::ModulePtr turbine(new Turbine("turbine", tjet.nodes_[3].get(),  tjet.nodes_[4].get(), tjet.mech_[0].get()));
    turbine->N2().Name_ = "Turbine Exit";					// Accesses the second node of the 'combustion' and sets its name.
    Turbine *turbineC = (Turbine *) turbine.get();				// TurbineC is a pointer to the Turbine class which is stored inside the pointer to the 'turbine' object - allows directed access of Turbine class.
    turbineC->efft_ = 0.85;
    turbineC->checkdesignpoint_ = ondesign;
	if (!ondesign){
		turbineC->Aratio_ = Aratio_;					
	}
	tjet.add(turbine);


    systemModule::ModulePtr nozzle(new AdaptedNozzle<EffIsenDuct>("nozzle", tjet.nodes_[4].get(), &tjet.N2()));
    AdaptedNozzle<EffIsenDuct> *nozzleC = (AdaptedNozzle<EffIsenDuct> *) nozzle.get();	// NozzleC is a pointer to the AdaptedNozzle template class which is stored inside the pointer to the 'nozzle' object - allows directed access of AdaptedNozzle template class.
  	nozzleC->etap_ = 1.0;
  	nozzleC->chokingMach_ = 1.0;			// Could this be changed - only works for up to M = 1.0
  	nozzleC->freeStream_ =  &inlet->N1();		// Sets the freestream conditions for the nozzle
  	tjet.add(nozzle);


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  																		Post Processing
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  	tjet.calculate();
  	//std::cout << freeStream.NfrX() << "\n";
	returned.sys = new systemModule(tjet);
	returned.thrust = tjet.thrust();


  	if (!((CaptureArea == false) && (ondesign == 0))) {

  		tjet.printOut(std::cout);
  	  	std::cout << "\n";
		fmt::print(std::cout, "{:<12} {:<12.2f}\n", "Thrust:", tjet.thrust());
		fmt::print(std::cout, "{:<12} {:<12.6f}\n", "kH:", turbineC->kH_);
		fmt::print(std::cout, "{:<12} {:<12.6f}\n", "Area Ratio:", turbineC->Aratio_);

  	}

    returned.T04 = tjet.nodes_[3]->getT0();
    returned.Aratio = turbineC->Aratio_;
    returned.kH = turbineC->kH_;

    returned.p9  = tjet.N2().getPress();
    returned.p05 = tjet.nodes_[4]->getp0();
    returned.p04 = tjet.nodes_[3]->getp0();
    returned.p03 = tjet.nodes_[1]->getp0();
    returned.p02 = tjet.nodes_[0]->getp0();
    returned.T02 = tjet.nodes_[0]->getT0();

    returned.V9  = tjet.N2().getU();
    returned.cp2 = tjet.nodes_[0]->Cp();
    returned.cp4 = tjet.nodes_[3]->Cp();
    returned.y5  = tjet.nodes_[4]->gamma();
    returned.y9  = tjet.N2().gamma();
    returned.f	 = (tjet.nodes_[2]->mfr() - tjet.nodes_[1]->mfr())/tjet.nodes_[1]->mfr();

	returned.min = tjet.nodes_[1]->mfr();
    returned.throttle = throttle;


    return {true, returned};


}
