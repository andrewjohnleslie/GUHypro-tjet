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
#include <backward/strstream>		// The header file 'strstream has been deprecated and replaced by 'stringstream'


#include "gasTurbine/Compressor.h"
#include "gasTurbine/Turbine.h"
#include "injectors/InjectionPhiNoTemp.h"
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
pyCycleTjetComparison(const double &p, const double &T, const double &Mach, const double &throttle, const int &verbosity, int &ondesign, double &T04_, double &Aratio_, double &kH_, double &CaptureArea) {

    outputData returned;

//	Node &freeStream =*new Node("/home/mark/git/GUHypro/test/pyCycleTjetMatching/nasa.xml", "mark");
//	Node &freeStream =*new Node("/home/mark/git/GUHypro/test/pyCycleTjetMatching/airNASA9.xml", "airNASA9");
//	Node &freeStream =*new Node("/home/mark/git/GUHypro/test/pyCycleTjetMatching/gri30_highT_modified_jetA.xml", "gri30");

	Node &freeStream =*new Node("/home/andrew/GUHypro-tjet/test/pyCycleTjetMatching/airNASA9_modified.yaml", "airNASA9");
//        Node &freeStream =*new Node("/home/andrew/GUHypro-tjet/test/TurboJet/gri30_highT_modified_jetA.yaml", "gri30");	

    freeStream.Name_ = "Freestream";

    Node &exit = *new Node(freeStream);
    exit.Name_ = "Nozzle exit to FS";

    Cantera::compositionMap initialComposition;
    initialComposition.emplace("O2", 0.209476);
    initialComposition.emplace("N2", 0.78084);
    initialComposition.emplace("Ar", 0.009365);
    initialComposition.emplace("CO2", 0.000319);


    freeStream.setTPX(T, p, initialComposition);
    vector<double> empty_thermos(freeStream.X().size());
//    freeStream.setU(Mach*freeStream.geta());
    freeStream.setU(0.0003404616);

    if ((ondesign == 0) && (Mach > 0.999)){
       	freeStream.setU(Mach*freeStream.geta());
   	}


    systemModule tjet("system", &freeStream, &exit);
    tjet.verbosity_ = verbosity;

    const int numNodes = 9;
    static vector<double> A(numNodes) ;

    A =  { 160551.799690, //16.0558,
    		0.32975864306204,
		    0.0933036424246,
    		0.14152630920488,
    		0.14152630920488,
    		0.14152630920488,
			0.14152630920488,
		    0.25472923957276,
		    0.18365164172824};



    // Set the areas of all the nodes
    const int numInternalNodes = numNodes - 2;

    if (CaptureArea){
//         freeStream.A(A[0]*CaptureArea);
//         freeStream.A(155265.95);
		 freeStream.A(0.198175);

     }
     else{
         freeStream.A(A[0]);
     }
    tjet.initNodes(numInternalNodes, freeStream);
    for (int i = 0; i < numInternalNodes; i++) {
    	tjet.nodes_[i]->A(A[i + 1]);
    }
    exit.A(A[numNodes - 1]);


    const int numShafts = 1;
    tjet.initShaft(numShafts);


/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Set up of Engine Configuration
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/

//    New here
    systemModule::ModulePtr inlet(new EffIsenDuct("inlet", &tjet.N1(), tjet.nodes_[0].get()));
    inlet->N2().Name_ = "Compressor Inlet";
    EffIsenDuct *inletC = (EffIsenDuct *) inlet.get();
    inletC->etap_ = 1.0;
	inletC->etaT_ = 1.0;
	inletC->choked_ = false;
	if ((Mach > 0.999) && (ondesign == 0)){
		inletC->choked_ = true;
	}
	tjet.add(inlet);

//    systemModule::ModulePtr inlet(new AdaptedInlet("inlet", &tjet.N1(), tjet.nodes_[0].get()));
//	AdaptedInlet *inletC = (AdaptedInlet *) inlet.get();
//	inlet->N2().Name_ = "LP Compressor Inlet";
//    inletC->chokingFeedback_ = inlet;
//	tjet.add(inlet);


    systemModule::ModulePtr compressor(new Compressor("compressor", tjet.nodes_[0].get(), tjet.nodes_[1].get(), tjet.mech_[0].get()));
    compressor->N2().Name_ = "Compressor Exit";
    Compressor *compressorC = (Compressor *) compressor.get();
	tjet.mech_[0]->eff_mech_ = 1.0;
    compressorC->effc_ = 0.87830073;
    compressorC->pi_c_ = 13.500000000;
//    compressorC->pi_c_ = 13.4995679;
    compressorC->checkdesignpoint_ = ondesign;
    if (ondesign == 0){
			compressorC->setT04(T04_);
			compressorC->kH_ = kH_;
			tjet.mech_[0]->Cp2_ =  1167.5;
			tjet.mech_[0]->Cp_ =  1226.51;
			tjet.mech_[0]->Ht_ =   403726.7*freeStream.mfr();
			tjet.mech_[0]->mfr_ =  freeStream.mfr();
			tjet.mech_[0]->kH_ = kH_;
//			compressorC->effc_ = 0.880705;
			compressorC->effc_ = 0.870492;

		}
	tjet.add(compressor);


    systemModule::ModulePtr injection(new InjectionPhiNoTemp<balanceMach>("injection", tjet.nodes_[1].get(),  tjet.nodes_[2].get()));
    injection->N2().Name_ = "Post Injection";
    InjectionPhiNoTemp<balanceMach> *injectionC = (InjectionPhiNoTemp<balanceMach> *) injection.get();
    injectionC->Phi_ = 1.0 * throttle;
    injectionC->PhiMax_ = 1.0;
    injectionC->chokingFeedback_ = injection;
    injectionC->chokingMach_ = 0.98;
    vector<double> Xf = empty_thermos;
    Xf[5] = 1;
	injectionC->Rcoeff_[0] = 71;
	injectionC->Rcoeff_[5] = 4;
	injectionC->Xf_ = Xf;
    tjet.add(injection);

    systemModule::ModulePtr duct1(new EffIsenDuct("Duct1", tjet.nodes_[2].get(),  tjet.nodes_[3].get()));
    duct1->N2().Name_ = "Combustion Inlet";
    EffIsenDuct *duct1C = (EffIsenDuct *) duct1.get();
    duct1C->etap_ = 1.0;
    duct1C->etaT_ = 1.0;
    duct1C->choked_ = false;
   	tjet.add(duct1);

    systemModule::ModulePtr combustion(new EffComb<balanceMach>("combustor", tjet.nodes_[3].get(),  tjet.nodes_[4].get()));
    combustion->N2().Name_ = "Combustion End";
    EffComb<balanceMach> *combustionC = (EffComb<balanceMach> *) combustion.get();
    combustionC->Rcoeff_[5] = 4;
	combustionC->Rcoeff_[0]  = 71;
	combustionC->Pcoeff_[2] = 48;
	combustionC->Pcoeff_[1]  = 46;
	combustionC->eta_ = 1;
	combustionC->Xswitch_ = 0;
	combustionC->chokingFeedback_ = injection;
	tjet.add(combustion);



    systemModule::ModulePtr duct2(new EffIsenDuct("Duct2", tjet.nodes_[4].get(),  tjet.nodes_[5].get()));
	duct2->N2().Name_ = "Turbine Inlet";
	EffIsenDuct *duct2C = (EffIsenDuct *) duct2.get();
	duct2C->etap_ = 0.982884;
	duct2C->etaT_ = 1;
	duct2C->choked_ = false;
	tjet.add(duct2);


    systemModule::ModulePtr turbine(new Turbine("turbine", tjet.nodes_[5].get(),  tjet.nodes_[6].get(), tjet.mech_[0].get()));
    turbine->N2().Name_ = "Turbine Exit";
    Turbine *turbineC = (Turbine *) turbine.get();
    turbineC->efft_ = 0.83933;//9031;
    turbineC->Xswitch_ = 0;
    turbineC->checkdesignpoint_ = ondesign;
	if (ondesign == 0){
    		turbineC->Aratio_ = Aratio_;
    	}
	tjet.add(turbine);


    systemModule::ModulePtr nozzle(new AdaptedNozzle<EffIsenDuct>("nozzle", tjet.nodes_[6].get(), &tjet.N2()));
    AdaptedNozzle<EffIsenDuct> *nozzleC = (AdaptedNozzle<EffIsenDuct> *) nozzle.get();
  	nozzleC->etap_ = 1.0;
  	nozzleC->chokingMach_ = 0.99;
  	nozzleC->freeStream_ =  &inlet->N1();
  	tjet.add(nozzle);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  	tjet.calculate();
  	tjet.printOut(std::cout);


  	std::cout << std::endl << "Thrust is: " << tjet.thrust() << " N" << std::endl;
    returned.sys = new systemModule(tjet);
    returned.thrust = tjet.thrust();
    returned.T04 = tjet.nodes_[5]->getT0();
    returned.Aratio = turbineC->Aratio_;
    returned.kH = turbineC->kH_;
    returned.throttle = throttle;

    returned.p9  = tjet.N2().getPress();
    returned.p05 = tjet.nodes_[6]->getp0();
    returned.p04 = tjet.nodes_[5]->getp0();
    returned.p03 = tjet.nodes_[1]->getp0();
    returned.p02 = tjet.nodes_[0]->getp0();
    returned.T02 = tjet.nodes_[0]->getT0();

    returned.V9  = tjet.N2().getU();

    returned.cp2 = tjet.nodes_[0]->Cp();
    returned.cp4 = tjet.nodes_[5]->Cp();
    returned.y5  = tjet.nodes_[6]->gamma();
    returned.y9  = tjet.N2().gamma();
    returned.f	 = (tjet.nodes_[2]->mfr() - tjet.nodes_[1]->mfr())/tjet.nodes_[1]->mfr();
//    returned.Ufs = tjet.N1().getU();


//    if (ondesign == 1) {
    		returned.min = tjet.nodes_[1]->mfr();
//    }




    return {true, returned};


}
