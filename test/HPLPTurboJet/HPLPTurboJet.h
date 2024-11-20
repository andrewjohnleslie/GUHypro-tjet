//
// Created by Mark De Luca on 14/12/18.
//
// For more information, contact m.de-luca.1@research.gla.ac.uk
//

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
HPLPTurboJet(const double &p, const double &T, const double &Mach, const double &throttle, const int &verbosity, int &ondesign, double &CaptureArea, double &T04_, double &T045_, double &Aratio_HP,double &Aratio_LP, double &kH_HP, double &kH_LP) {

    outputData returned;

//    Node &freeStream =*new Node("/home/andrew/GUHypro-tjet/test/TurboJet/gri30_highT_modified_jetA.yaml", "gri30");
    Node &freeStream =*new Node("/home/andrew/GUHypro-tjet/test/pyCycleTjetMatching/airNASA9.yaml", "airNASA9");

    freeStream.Name_ = "Freestream";
    
    Node &exit = *new Node(freeStream);
    exit.Name_ = "Nozzle Exit";

    Cantera::compositionMap initialComposition;
    initialComposition.emplace("O2", 0.209);
    initialComposition.emplace("N2", 1 - 0.209);

    freeStream.setTPX(T, p, initialComposition);
    vector<double> empty_thermos(freeStream.X().size());
    freeStream.setU(Mach * freeStream.geta());

    systemModule twospooltjet("twospooltjet", &freeStream, &exit);
    twospooltjet.verbosity_ = verbosity;

    const int numNodes = 23;
//	const double AA = 1.0;
	static vector<double> A(numNodes);
//    A = {AA, AA, AA, AA, AA, AA, AA, AA, AA*0.9};
//    A = {1.0946, 0.3690, 0.3620, 0.2366, 0.3194, 0.4352, 0.4136, 0.8970, 1.0399, 0.7727, 1.3908};		// An area ratio that assumes a denominator of ~1550 in^2 


//    A = {1.380535, 1.0946, 0.36895, 0.3620, 0.2366, 0.2106, 0.3194, 0.3194, 0.3194, 0.3194, 0.4352, 0.4352, 0.4352, 0.41364, 0.89704, 1.0399, 0.0260, 0.0236,  0.0118,  0.0118, 0.0024, 0.77265, 1.3908};
    A = {1.380535, 1.0946, 0.36895, 0.3620, 0.2366, 0.2106, 0.3194, 0.3194, 0.3194, 0.3194, 0.4352, 0.4352, 0.4352, 0.41364, 0.89704, 1.0399, 0.0260, 0.0236,  0.0118,  0.0118, 0.0024, 0.773, 1.3908};
//          N1       0       1        2	       3       4       5    	6       7		8	  	 9      10	    11		12		 13		 14       15	 16		  17		18	    19		20		N2
    // Set the areas of all the nodes

    static vector<string> Names = {
          	"LP Compressor Inlet",			//[0]
          	"LP Compressor Exit",			//[1]
   			"HP Compressor Inlet",			//[2]
   			"HP Compressor Exit",			//[3]
   			"Core Duct Inlet",				//[4]
   			"Injection Inlet",				//[5]
          	"Combustion Inlet",				//[6]
   			"Combustion Exit",				//[7]
			"HP Turbine Inlet",				//[8]
			"NGV Cooling",					//[9]
   			"HP Turbine Inter",				//[10]
   			"HP Turbine Exit",				//[11]
   		    "LP Turbine Inlet",				//[12]
   			"LP Turbine Exit",				//[13]
   			"Gas Generator Exit",			//[14]
   			"Bleed Air",					//[15]
			"Total Cooling Air",			//[16]
			"Cooling Air1",					//[17]
			"Cooling Air2",					//[18]
			"Cabin Bleed",					//[19]
			"Nozzle Throat"					//[20]
      	};


    const int numInternalNodes = numNodes - 2;
    freeStream.A(A[0]);

    if (CaptureArea){
        freeStream.A(A[0]*CaptureArea);
    }
    else{
        freeStream.A(A[0]);
    }

    twospooltjet.initNodes(numInternalNodes, freeStream);
    for (int i = 0; i < numInternalNodes; i++) {
    	twospooltjet.nodes_[i]->A(A[i + 1]);
    	twospooltjet.nodes_[i]->Name_ = Names[i];

    }
    exit.A(A[numNodes - 1]);
//    exit.Amax_ = A[numNodes];


    twospooltjet.initShaft(2);

    /*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Set up of Engine Configuration
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/

//    if (ondesign == 0 && Mach <1.00){
		systemModule::ModulePtr inlet(new AdaptedInlet("inlet", &twospooltjet.N1(), twospooltjet.nodes_[0].get()));
		AdaptedInlet *inletC = (AdaptedInlet *) inlet.get();
		twospooltjet.add(inlet);
//    }
//    else {
//		systemModule::ModulePtr inlet(new Inlet("inlet", &twospooltjet.N1(), twospooltjet.nodes_[0].get()));
//		Inlet *inletC = (Inlet *) inlet.get();
//		twospooltjet.add(inlet);
//    }

    systemModule::ModulePtr lpcompressor(new Compressor("lpcompressor", twospooltjet.nodes_[0].get(), twospooltjet.nodes_[1].get(), twospooltjet.mech_[0].get()));
    Compressor *lpcompressorC = (Compressor *) lpcompressor.get();
    lpcompressorC->effc_ = 0.8782;
    lpcompressorC->pi_c_ = 4.1;
    lpcompressorC->checkdesignpoint_ = ondesign;
    twospooltjet.mech_[0]->eff_mech_ = 0.99;		// NOT IN TURBOJET.H
    if (ondesign == 0){
    	lpcompressorC->setT04(T045_);
    	lpcompressorC->kH_ = kH_LP;
//        lpcompressorC->effc_ = 0.885; //For Mach = 1.9;
//        lpcompressorC->effc_ = 0.827; //For Mach = 0.3241, Alt = 0;
//        lpcompressorC->effc_ = 0.758; //For Mach = 1.3,    Alt = 10668;
        twospooltjet.mech_[0]->Cp2_ =  twospooltjet.nodes_[0]->Cp();
    	twospooltjet.mech_[0]->Cp_ =  twospooltjet.nodes_[0]->Cp();
    	twospooltjet.mech_[0]->mfr_ =  twospooltjet.nodes_[0]->mfr();
    	twospooltjet.mech_[0]->kH_ = kH_LP;
    	twospooltjet.mech_[0]->Ht_ = 184286*freeStream.mfr();
		}
//    lpcompressorC->chokingFeedback_ = inlet;
    twospooltjet.add(lpcompressor);


	systemModule::ModulePtr ductlphpcomp(new EffIsenDuct("ductlphpcomp", twospooltjet.nodes_[1].get(), twospooltjet.nodes_[2].get()));
	EffIsenDuct *ductlphpcompC = (EffIsenDuct *) ductlphpcomp.get();
	ductlphpcompC->etap_ = 0.99;
	twospooltjet.add(ductlphpcomp);


    systemModule::ModulePtr hpcompressor(new Compressor("hpcompressor", twospooltjet.nodes_[2].get(), twospooltjet.nodes_[3].get(), twospooltjet.mech_[1].get()));
    Compressor *hpcompressorC = (Compressor *) hpcompressor.get();
    hpcompressorC->effc_ = 0.8402;
    hpcompressorC->pi_c_ = 2.9;
    hpcompressorC->checkdesignpoint_ = ondesign;
    twospooltjet.mech_[1]->eff_mech_ = 0.99;
    if (ondesign == 0){
    	hpcompressorC->setT04(T04_);
    	hpcompressorC->kH_ = kH_HP;
//        hpcompressorC->effc_ = 0.890;  //For Mach = 0.3241, Alt = 0;
//        hpcompressorC->effc_ = 0.9;    //For Mach = 1.3,    Alt = 10668;
    	twospooltjet.mech_[1]->Cp2_ =  twospooltjet.nodes_[0]->Cp();
    	twospooltjet.mech_[1]->Cp_ =  twospooltjet.nodes_[0]->Cp();
    	twospooltjet.mech_[1]->mfr_ =  twospooltjet.nodes_[0]->mfr();
    	twospooltjet.mech_[1]->kH_ = kH_HP;
    	twospooltjet.mech_[1]->Ht_ = 267345*freeStream.mfr();
		}
    twospooltjet.add(hpcompressor);


	systemModule::ModulePtr bleed1(new simpleSplitter("bleed1", twospooltjet.nodes_[3].get(), twospooltjet.nodes_[4].get(), twospooltjet.nodes_[15].get()));
	simpleSplitter *bleed1C = (simpleSplitter *) bleed1.get();
//	bleed1C->etap_ = 1.0;
	twospooltjet.add(bleed1);

	systemModule::ModulePtr split(new simpleSplitter("split", twospooltjet.nodes_[15].get(), twospooltjet.nodes_[16].get(), twospooltjet.nodes_[19].get()));
	simpleSplitter *splitC = (simpleSplitter *) split.get();
//	bleed1C->etap_ = 1.0;
	twospooltjet.add(split);

	systemModule::ModulePtr split2(new simpleSplitter("split2", twospooltjet.nodes_[16].get(), twospooltjet.nodes_[17].get(), twospooltjet.nodes_[18].get()));
	simpleSplitter *split2C = (simpleSplitter *) split2.get();
//	bleed1C->etap_ = 1.0;
	twospooltjet.add(split2);

	systemModule::ModulePtr duct2(new EffIsenDuct("duct2", twospooltjet.nodes_[4].get(), twospooltjet.nodes_[5].get()));
	EffIsenDuct *duct2C = (EffIsenDuct *) duct2.get();
	duct2C->etap_ = 0.98;
	twospooltjet.add(duct2);


    systemModule::ModulePtr injection(new InjectionPhiNoTemp<balanceMach>("injection", twospooltjet.nodes_[5].get(),  twospooltjet.nodes_[6].get()));
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
    twospooltjet.add(injection);


    systemModule::ModulePtr combustion(new EffComb<balanceMach>("combustor", twospooltjet.nodes_[6].get(),  twospooltjet.nodes_[7].get()));
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
	twospooltjet.add(combustion);


//	systemModule::ModulePtr duct2(new EffIsenDuct("duct2", twospooltjet.nodes_[5].get(), twospooltjet.nodes_[6].get()));
//	EffIsenDuct *duct2C = (EffIsenDuct *) duct2.get();
//	duct2->N2().Name_ = "HP Turbine Inlet";
//	duct2C->etap_ = 0.98;
//	twospooltjet.add(duct2);


	systemModule::ModulePtr duct5(new EffIsenDuct("duct5", twospooltjet.nodes_[7].get(), twospooltjet.nodes_[8].get()));
	EffIsenDuct *duct5C = (EffIsenDuct *) duct5.get();
	duct5C->etap_ = 0.99435;
	twospooltjet.add(duct5);


	systemModule::ModulePtr recombinebleed1(new Mixer("recombinebleed1", twospooltjet.nodes_[8].get(), twospooltjet.nodes_[9].get(), twospooltjet.nodes_[17].get()));
	Mixer *recombinebleed1C = (Mixer *) recombinebleed1.get();
	recombinebleed1C->verbosity_ = 0;
	recombinebleed1C->teta_ = 42.2;
//	recombinebleed1C->chokingFeedback_ = inlet;
//	recombinebleed1C->eff_ = 0.82;
	twospooltjet.add(recombinebleed1);


    systemModule::ModulePtr hpturbine(new Turbine("hpturbine", twospooltjet.nodes_[9].get(),  twospooltjet.nodes_[10].get(), twospooltjet.mech_[1].get()));
    Turbine *hpturbineC = (Turbine *) hpturbine.get();
    hpturbineC->efft_ = 0.8785;
//    hpturbineC->efft_ = 0.91325;
    hpturbineC->checkdesignpoint_ = ondesign;
	if (ondesign == 0){
    		hpturbineC->Aratio_ = Aratio_HP;
    		hpturbineC->kH_ = kH_HP;
    	    hpturbineC->efft_ = 0.8785;
	}
    twospooltjet.add(hpturbine);

	systemModule::ModulePtr recombinebleed2(new Mixer("recombinebleed2", twospooltjet.nodes_[10].get(), twospooltjet.nodes_[11].get(), twospooltjet.nodes_[18].get()));
	Mixer *recombinebleed2C = (Mixer *) recombinebleed2.get();
	recombinebleed2C->verbosity_ = 0;
	recombinebleed2C->teta_ = 75.55;
//	recombinebleed1C->chokingFeedback_ = inlet;
//	recombinebleed2C->eff_ = 0.82; // Nozzle Throat Area = 0.7756
	twospooltjet.add(recombinebleed2);


	systemModule::ModulePtr duct3(new EffIsenDuct("duct3", twospooltjet.nodes_[11].get(), twospooltjet.nodes_[12].get()));
	EffIsenDuct *duct3C = (EffIsenDuct *) duct3.get();
	duct3C->etap_ = 0.98;
	twospooltjet.add(duct3);


    systemModule::ModulePtr lpturbine(new Turbine("lpturbine", twospooltjet.nodes_[12].get(),  twospooltjet.nodes_[13].get(), twospooltjet.mech_[0].get()));
    Turbine *lpturbineC = (Turbine *) lpturbine.get();
//    lpturbineC->efft_ = 0.922; //Nozzle throat area is 0.7725
    lpturbineC->efft_ = 0.8890; //Nozzle throat area is 0.79981
    lpturbineC->checkdesignpoint_ = ondesign;
	if (ondesign == 0){
    		lpturbineC->Aratio_ = Aratio_LP;
    		lpturbineC->kH_ = kH_LP;
    	    lpturbineC->efft_ = 0.8890;
	}
    twospooltjet.add(lpturbine);


	systemModule::ModulePtr duct4(new EffIsenDuct("duct4", twospooltjet.nodes_[13].get(), twospooltjet.nodes_[14].get()));
	EffIsenDuct *duct4C = (EffIsenDuct *) duct4.get();
	duct4C->etap_ = 0.98;
	twospooltjet.add(duct4);

	systemModule::ModulePtr ConvergentNozzle(new EffIsenDuct("ConvNozzle", twospooltjet.nodes_[14].get(), twospooltjet.nodes_[20].get()));
	EffIsenDuct *ConvergentNozzleC = (EffIsenDuct *) ConvergentNozzle.get();
//	ConvergentNozzleC->verbosity_ = 1;
	ConvergentNozzleC->chokingMach_= 1.0;
	ConvergentNozzleC->choked_ = 1;
	twospooltjet.add(ConvergentNozzle);

	systemModule::ModulePtr DivNozzle(new EffIsenDuct("DivNozzle", twospooltjet.nodes_[20].get(),  &twospooltjet.N2()));
	EffIsenDuct *DivNozzleC = (EffIsenDuct *) DivNozzle.get();
//	ConvNozzleC->etap_ = 0.98;
	twospooltjet.add(DivNozzle);

//    systemModule::ModulePtr nozzle(new AdaptedNozzle<EffIsenDuct>("nozzle", twospooltjet.nodes_[14].get(), &twospooltjet.N2()));
//    AdaptedNozzle<EffIsenDuct> *nozzleC = (AdaptedNozzle<EffIsenDuct> *) nozzle.get();
//  	nozzleC->etap_ = 1.0;
//  	nozzleC->chokingMach_ = 0.99;
//  	nozzleC->freeStream_ =  &freeStream;
//  	twospooltjet.add(nozzle);




  	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  	//// Execution and Post Processing
  	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



  	twospooltjet.calculate();
  	returned.sys = new systemModule(twospooltjet);
    returned.thrust =twospooltjet.thrust();
    returned.Mfr =twospooltjet.propellantMfr();

  	if (!((CaptureArea == false) && (ondesign == 0))) {

  		twospooltjet.printOut(std::cout);
  	  	std::cout << "\n" << "Thrust is: " << twospooltjet.thrust() << " N" << "\n";

	}




    returned.T04 = combustion->N2().getT0();
   	returned.Aratio_HP = hpturbineC->Aratio_;
   	returned.kH_HP = hpturbineC->kH_;
    returned.T045 = hpturbine->N2().getT0();
    returned.Aratio_LP = lpturbineC->Aratio_;
    returned.kH_LP = lpturbineC->kH_;

    returned.p9  = twospooltjet.N2().getPress();
    returned.p05 = lpturbine->N2().getp0();
    returned.p04 = combustion->N2().getp0();
    returned.p03 = hpcompressor->N2().getp0();
    returned.p02 = lpcompressor->N1().getp0();
    returned.T02 = lpcompressor->N1().getT0();

    returned.Ufs = twospooltjet.N1().getU();
    returned.V9  = twospooltjet.N2().getU();

    returned.cp2 = lpcompressor->N1().Cp();
    returned.cp4 = combustion->N2().Cp();
    returned.y5  = lpturbine->N2().gamma();
    returned.y9  = twospooltjet.N2().gamma();
    returned.f	 = (injectionC->N2().mfr() - injectionC->N1().mfr())/injectionC->N1().mfr();

    returned.Mach = Mach;
    returned.throttle = throttle;

	returned.min = twospooltjet.nodes_[1]->mfr();

    //  //Prints Chemical Composition of each Node
	//	std::cout << twospooltjet.N1().Name_ << "\n";
	//	std::cout << twospooltjet.N1().report() << "\n";
	//
	//	for (std::size_t i = 0; i != numInternalNodes; i++){
	//	    std::cout << "i points to: " << twospooltjet.nodes_[i]->Name_ << "\n";
	//		std::cout << twospooltjet.nodes_[i]->report() << "\n";;
	//
	//	}
	//
	//	std::cout << twospooltjet.N2().Name_ << "\n";
	//	std::cout << twospooltjet.N2().report() << "\n";

    return {true, returned};


}
