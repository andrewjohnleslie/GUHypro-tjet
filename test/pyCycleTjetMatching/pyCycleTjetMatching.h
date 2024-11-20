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
pyCycleTjetMatching(const double &p, const double &T, const double &Mach, const double &throttle, const int &verbosity, int &ondesign, double &T04_, double &Aratio_, double &kH_, double &CaptureArea) {

    outputData returned;

//	Node &freeStream =*new Node("/home/mark/git/GUHypro/test/pyCycleTjetMatching/nasa.xml", "mark");
//	Node &freeStream =*new Node("/home/mark/git/GUHypro/test/pyCycleTjetMatching/airNASA9.xml", "airNASA9");
	Node &freeStream =*new Node("/home/andrew/GUHypro-tjet/test/pyCycleTjetMatching/airNASA9_modified.yaml", "airNASA9");
//      Node &freeStream =*new Node("/home/andrew/GUHypro-tjet/test/TurboJet/gri30_highT_modified_jetA.yaml", "gri30");

    freeStream.Name_ = "Freestream";

    Node &exit = *new Node(freeStream);
    exit.Name_ = "Nozzle exit to FS";

    Cantera::compositionMap initialComposition;
    initialComposition.emplace("O2", 0.209476);
    initialComposition.emplace("N2", 0.78084);
    initialComposition.emplace("Ar", 0.009365);
    initialComposition.emplace("CO2", 0.000319);

    std::vector<double> pyCycleTurbInXDP (freeStream.X().size());
    std::vector<double> pyCycleTurbOutXDP (freeStream.X().size());
    std::vector<double> pyCycleTurbInXOD0 (freeStream.X().size());
	std::vector<double> pyCycleTurbOutXOD0 (freeStream.X().size());

    pyCycleTurbInXDP  = {0.152134831117542, 0.0346687377373258, 0.0364914697290112, 0.767116912136294, 0.00920270079087897, 0, 0.000370007095074531, 1.16187917087244e-05, 3.6293843787906e-06};
    pyCycleTurbOutXDP = {0.152315752009611, 0.0346704142289617, 0.036491334525617, 0.767289899340093, 0.00920266602059696, 0, 2.7588240861387e-05, 2.26350368790788e-06, 5.02663151113896e-08};
//    pyCycleTurbInXOD0  = {0.155059919047964, 0.0329207837841992, 0.0346670726758753, 0.767841493657197, 0.0092108867503504, 0, 0.000287056619261036, 1.04726988577057e-05, 2.25018327019577e-06};
//    pyCycleTurbOutXOD0 = {0.155202342523783, 0.032921785595832, 0.0346669444287499, 0.767976817335059, 0.00921085267566059, 0, 1.92365457258093e-05, 1.96322112491362e-06, 2.58099407924636e-08};
    pyCycleTurbInXOD0  = {0.165554010379013, 0.0265884748945771, 0.0280596762701859, 0.770369617864911, 0.00924054417450375, 0, 0.000176520356100815, 1.03461092116786e-05 , 7.6893990123626e-07 };
	pyCycleTurbOutXOD0 = {0.1656456649442, 0.0265887521639927, 0.0280595614017686, 0.770454316654423, 0.00924050634631663, 0, 9.37910035050147e-06, 1.78166816718895e-06, 5.85712015083846e-09};

    freeStream.setTPX(T, p, initialComposition);
    vector<double> empty_thermos(freeStream.X().size());
    freeStream.setU(0.0003404616);

    if ((ondesign == 0) && (Mach > 0.999)){
    	freeStream.setU(Mach*freeStream.geta());
	}


    systemModule tjet("system", &freeStream, &exit);
    tjet.verbosity_ = verbosity;

    const int numNodes = 9;
    static vector<double> A(numNodes) ;

    A =  {160551.79969,
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

    if (ondesign == 0){
//            freeStream.A(155265.95 );
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
 
	std::cout << "TESTING - pycycletjetmatch.h - LINE 147 " << std::endl;
    systemModule::ModulePtr inlet(new EffIsenDuct("inlet", &tjet.N1(), tjet.nodes_[0].get()));
    inlet->N2().Name_ = "Compressor Inlet";
    EffIsenDuct *inletC = (EffIsenDuct *) inlet.get();
    inletC->etap_ = 1.0;
	inletC->etaT_ = 1.0;
	inletC->choked_ = false;
    if ((ondesign == 0) && (Mach > 0.999)){
    	inletC->choked_ = true;

    }

	tjet.add(inlet);


    systemModule::ModulePtr compressor(new Compressor("compressor", tjet.nodes_[0].get(), tjet.nodes_[1].get(), tjet.mech_[0].get()));
    compressor->N2().Name_ = "Compressor Exit";
    Compressor *compressorC = (Compressor *) compressor.get();
	tjet.mech_[0]->eff_mech_ = 1.0;
    compressorC->effc_ = 0.87830073;
    compressorC->pi_c_ = 13.500000000;
    compressorC->checkdesignpoint_ = ondesign;
    if (ondesign == 0){
			compressorC->setT04(T04_);
			compressorC->kH_ = kH_;
			tjet.mech_[0]->Cp2_ =  1228;
			tjet.mech_[0]->Cp_ =  1228;
			tjet.mech_[0]->Ht_ =  24000000;
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


//    systemModule::ModulePtr injection(new InjectionPhi<balanceMach>("injection", tjet.nodes_[1].get(),  tjet.nodes_[2].get()));
//    injection->N2().Name_ = "Post Injection";
//    InjectionPhi<balanceMach> *injectionC = (InjectionPhi<balanceMach> *) injection.get();
//    injectionC->Phi_ = 1.0 * throttle;
//    injectionC->PhiMax_ = 1.0;
//    injectionC->chokingFeedback_ = injection;
//    injectionC->chokingMach_ = 0.98;
//    vector<double> Xf = empty_thermos;
//    injectionC->Tf_ = 841.37;
//    if (ondesign == 0){
//    	injectionC->Tf_ = 846;
//    }
//    Xf[5] = 1;
//	injectionC->Rcoeff_[0] = 71;
//	injectionC->Rcoeff_[5] = 4;
//	injectionC->Xf_ = Xf;
//    tjet.add(injection);

	std::cout << "TESTING - pycycletjetmatch.h - LINE 217 " << std::endl;
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
	combustionC->Xswitch_ = 1;
	combustionC->XSET_ = pyCycleTurbInXDP;
	if (ondesign == 0){
		combustionC->XSET_ = pyCycleTurbInXOD0;
	}
	combustionC->chokingFeedback_ = injection;
	tjet.add(combustion);

	std::cout << "TESTING - pycycletjetmatch.h - LINE 242 " << std::endl;

    systemModule::ModulePtr duct2(new EffIsenDuct("Duct2", tjet.nodes_[4].get(),  tjet.nodes_[5].get()));
	duct2->N2().Name_ = "Turbine Inlet";
	EffIsenDuct *duct2C = (EffIsenDuct *) duct2.get();
	duct2C->etap_ = 0.982884;
	if (ondesign == 0){
		duct2C->etap_ = 0.9828;
	}
	duct2C->etaT_ = 1;
	duct2C->choked_ = false;
	tjet.add(duct2);


    systemModule::ModulePtr turbine(new Turbine("turbine", tjet.nodes_[5].get(),  tjet.nodes_[6].get(), tjet.mech_[0].get()));
    turbine->N2().Name_ = "Turbine Exit";
    Turbine *turbineC = (Turbine *) turbine.get();
    turbineC->efft_ = 0.838845;
    turbineC->Xswitch_ = 1;
    turbineC->XSET_ = pyCycleTurbOutXDP;
    turbineC->checkdesignpoint_ = ondesign;
	if (ondesign == 0){
    		turbineC->Aratio_ = Aratio_;
    	    turbineC->XSET_ = pyCycleTurbOutXOD0;
//    	    turbineC->efft_ = 0.838;
	}
	tjet.add(turbine);

	std::cout << "TESTING - pycycletjetmatch.h - LINE 270 " << std::endl;

    systemModule::ModulePtr nozzle(new AdaptedNozzle<EffIsenDuct>("nozzle", tjet.nodes_[6].get(), &tjet.N2()));
    AdaptedNozzle<EffIsenDuct> *nozzleC = (AdaptedNozzle<EffIsenDuct> *) nozzle.get();
  	nozzleC->etap_ = 1.0;
  	nozzleC->chokingMach_ = 0.99;
  	nozzleC->freeStream_ =  &inlet->N1();
  	tjet.add(nozzle);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	std::cout << "TESTING - pycycletjetmatch.h - LINE 281 " << std::endl;
  	tjet.calculate();	// Code fails here!!
  	std::cout << "TESTING - pycycletjetmatch.h - LINE 283 " << std::endl;
	tjet.printOut(std::cout);

	std::cout << "TESTING - pycycletjetmatch.h - LINE 286 " << std::endl;

// 	if (!((CaptureArea == false) && (ondesign == 0))) {

	std::cout << "\n" << "Thrust is: " << tjet.thrust() << " N" << "\n";
	returned.min = tjet.nodes_[0]->mfr();

//	}

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

    if (ondesign == 1) {
    		returned.min = tjet.nodes_[1]->mfr();
    }




    return {true, returned};


}
