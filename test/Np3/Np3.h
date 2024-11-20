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
#include "gasTurbine/TurbineBleed.h"
#include "injectors/InjectionPhi.h"
#include "combustion/Combustion.h"
#include "combustion/EffComb.h"
#include "combustion/EffPhiComb.h"
#include "combustion/EquilCombustion.h"
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
#include "simpleSplitter.h"
#include "Mixer.h"
#include "ParMixer.h"
#include "mixerModule.h"
#include "NeutralLink.h"


#include "cantera/thermo/ThermoPhase.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/Phase.h"
#include <memory>
#include "typedefs.h"
#include "offDesignData.h"


using namespace hypro;
std::pair<bool, outputData>
Np3(const double &p, const double &T, const double &Mach, const double &throttle, const int &verbosity) {


    outputData returned;

    //Node &freeStream =*new Node("/home/mark/git/GUHypro/test/TurboJet/gri30_highT_modified_jetA.xml", "gri30");
    Node &freeStream =*new Node("/home/andrew/GUHypro-tjet/test/TurboJet/gri30_highT_modified_jetA.yaml", "gri30");
    freeStream.Name_ = "Freestream";


    Cantera::compositionMap initialComposition;
    initialComposition.emplace("O2", 0.209);
    initialComposition.emplace("N2", 1 - 0.209);

//    initialComposition.emplace("O2", 0.209476);
//    initialComposition.emplace("N2", 0.78084);
//    initialComposition.emplace("AR", 0.009365);
//    initialComposition.emplace("CO2",0.000319);
    // These values taken from pyCycle N+3 Data

    freeStream.setTPX(T, p, initialComposition);
    vector<double> empty_thermos(freeStream.X().size());
    freeStream.setU(Mach * freeStream.geta());

    Node &exit = *new Node(freeStream);
    exit.Name_ = "mixer to exit";// N2

    systemModule np3("system", &freeStream, &exit);
    np3.verbosity_ = verbosity;

    static vector<double> Areas = {
    		4.0961, 				// N1
			4.5870, 				//[0]
			4.5793, 				//[1]
			4.3961,			 		//[2]
			0.18329, 				//[3]
			0.1851,			 		//[4]
			0.073419,  				//[5]
			0.07458,  				//[6]
			0.07458*0.98,			//[7]
			0.07458*0.02,			//[8]
			0.011097,  				//[9]
			0.011097*0.8653,		//[10]
			0.011097*0.1347,		//[11]
			0.0096129, 				//[12]
			0.043548, 				//[13]
			0.043548, 				//[14]
			0.043548, 				//[15]
			0.059935, 				//[16]
			0.042774, 				//[17]
			0.44619, 				//[18]
			0.60968, 				//[19]
			0.2531, 				//[20]
			4.4630, 				//[21]
			3.0808, 				//[22]
			4.0951					//[N2]
    };


    static vector<string> Names = {
    		"fan entrance",			//[0]
    		"fan exit",				//[1]
    		"fan bypass",			//[2]
    		"fan core",				//[3]
    		"LPC entrance",			//[4]
    		"LPC exit",				//[5]
    		"LPC bleed entrance",	//[6]
    		"HPC entrance",			//[7]
    		"LPC bleed flow",		//[8]
    		"HPC bleed entrance",	//[9]
    		"HPC exit",				//[10]
    		"HPC bleed flow",		//[11]
    		"injection end",		//[12]
			"burner entrance",		//[13]
			"burner exit (inc)",	//[14]
			"burner exit (dp)",		//[15]
			"HPT exit",				//[16]
    		"LPT entrance",			//[17]
    		"LPT exit",				//[18]
    		"core nozzle entrance",	//[19]
    		"core nozzle exit",		//[20]
    		"fan nozzle entrance",	//[21]
    		"fan nozzle exit"		//[22]
	};

    static int numNodes = Areas.size();
    static int numInternalNodes = numNodes - 2;

    freeStream.A(Areas[0]);

	np3.initNodes(numInternalNodes, freeStream);
	for (int i = 0; i < numInternalNodes; i++) {
		np3.nodes_[i]->A(Areas[i + 1]);
		np3.nodes_[i]->Name_ = Names[i];
	}
	exit.A(Areas[numNodes - 1]);

    const int numShafts = 2;
	np3.initShaft(numShafts);
	np3.mech_[0]->eff_mech_ = 0.99;
	np3.mech_[1]->eff_mech_ = 0.9655;

	/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Set up of Engine Configuration
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/

    systemModule::ModulePtr inlet(new EffIsenDuct("inlet", &np3.N1(), np3.nodes_[0].get()));
    EffIsenDuct *inletC = (EffIsenDuct *) inlet.get();
    inletC->etap_ = 0.998;
	inletC->etaT_ = 1.0;
	inletC->choked_ = false;
	np3.add(inlet);

	systemModule::ModulePtr fan(new Compressor("fan", np3.nodes_[0].get(), np3.nodes_[1].get(), np3.mech_[0].get()));
	Compressor *fanC = (Compressor *) fan.get();
	fanC->effc_ = 0.97;
	fanC->pi_c_ = 1.3;
	np3.add(fan);

	systemModule::ModulePtr split(new simpleSplitter("splitter", np3.nodes_[1].get(), np3.nodes_[2].get(), np3.nodes_[3].get()));
	np3.add(split);


	systemModule::ModulePtr link1(new EffIsenDuct("linkfan2lpc", np3.nodes_[3].get(), np3.nodes_[4].get()));
	EffIsenDuct *link1C = (EffIsenDuct *) link1.get();
	link1C->etap_ = 0.99;
	link1C->etaT_ = 1.0;
	np3.add(link1);


	systemModule::ModulePtr lpcompressor(new Compressor("lpcompressor", np3.nodes_[4].get(), np3.nodes_[5].get(), np3.mech_[0].get()));
	Compressor *lpcompressorC = (Compressor *) lpcompressor.get();
	lpcompressorC->effc_ = 0.9050;
	lpcompressorC->pi_c_ = 3.0;
	np3.add(lpcompressor);


	systemModule::ModulePtr link2(new EffIsenDuct("linklpc2hpc", np3.nodes_[5].get(), np3.nodes_[6].get()));
	EffIsenDuct *link2C = (EffIsenDuct *) link2.get();
	link2C->etap_ = 0.985;
	link2C->etaT_ = 1.0;
	np3.add(link2);

	systemModule::ModulePtr bleed1(new simpleSplitter("bleed1", np3.nodes_[6].get(), np3.nodes_[7].get(), np3.nodes_[8].get()));
	np3.add(bleed1);

	systemModule::ModulePtr hpcompressor(new Compressor("hpcompressor", np3.nodes_[7].get(), np3.nodes_[9].get(), np3.mech_[1].get()));
	Compressor *hpcompressorC = (Compressor *) hpcompressor.get();
	hpcompressorC->effc_ = 0.89;
	hpcompressorC->pi_c_ = 14.103;
	np3.add(hpcompressor);

	systemModule::ModulePtr bleed2(new simpleSplitter("bleed2", np3.nodes_[9].get(), np3.nodes_[10].get(), np3.nodes_[11].get()));
	np3.add(bleed2);

	systemModule::ModulePtr injection(new InjectionPhi<balanceMach>("injection", np3.nodes_[10].get(),  np3.nodes_[12].get()));
	InjectionPhi<balanceMach> *injectionC = (InjectionPhi<balanceMach> *) injection.get();
	injectionC->Phi_ = 1.0 * throttle;
	injectionC->PhiMax_ = 1.0;
	injectionC->chokingFeedback_ = injection;
	injectionC->chokingMach_ = 0.98;
	injectionC->Tf_ = 836.41;
	vector<double> Xf = empty_thermos;
	Xf[53] = 1;
	injectionC->Rcoeff_[3] = 71;
	injectionC->Rcoeff_[53] = 4;
	injectionC->Xf_ = Xf;
//	injectionC->Cf_ = 0.0035;
//	injectionC->L_ = 0.1;
//	injectionC->D_ = 2 * sqrt(injection->N1().A()/(3.141));
//	injectionC->solver_verbosity = true;
	np3.add(injection);


	systemModule::ModulePtr diffuser(new EffIsenDuct("diffuser", np3.nodes_[12].get(), np3.nodes_[13].get()));
	EffIsenDuct *diffuserC = (EffIsenDuct *) diffuser.get();
	diffuserC->etap_ = 1.0;
	diffuserC->etaT_ = 1.0;
	np3.add(diffuser);



	systemModule::ModulePtr combustion(new EffComb<Friction>("combustor", np3.nodes_[13].get(),  np3.nodes_[14].get()));
	EffComb<Friction> *combustionC = (EffComb<Friction> *) combustion.get();
	combustionC->Rcoeff_[53] = 4;
	combustionC->Rcoeff_[3]  = 71;
	combustionC->Pcoeff_[15] = 48;
	combustionC->Pcoeff_[5]  = 46;
	combustionC->eta_ = 1.0;
	combustionC->Cf_ = 0;
	combustionC->D_ = 2 * sqrt(combustion->N1().A()/(3.141));
	combustionC->L_ = 1;
	combustionC->solver_verbosity = false;
	combustionC->chokingFeedback_ = injection;
	np3.add(combustion);

	systemModule::ModulePtr combustpdrop(new EffIsenDuct("combustpdrop", np3.nodes_[14].get(), np3.nodes_[15].get()));
	EffIsenDuct *combustpdropC = (EffIsenDuct *) combustpdrop.get();
	combustpdropC->etap_ = 0.965;
	combustpdropC->etaT_ = 1.0;
	np3.add(combustpdrop);

//	systemModule::ModulePtr recombinebleed1(new Mixer("recombinebleed1", np3.nodes_[15].get(), np3.nodes_[16].get(), np3.nodes_[11].get()));
//	Mixer *recombinebleed1C = (Mixer *) recombinebleed1.get();
//	recombinebleed1C->verbosity_ = 0;
////	recombinebleed1C->chokingFeedback_ = inlet;
//	recombinebleed1C->eff_ = 0.99;
//	np3.add(recombinebleed1);



    systemModule::ModulePtr hpturbine(new TurbineBleed("hpturbine", np3.nodes_[15].get(),  np3.nodes_[16].get(), np3.nodes_[11].get(), np3.mech_[1].get()));
    TurbineBleed *hpturbineC = (TurbineBleed *) hpturbine.get();
    hpturbineC->efft_ = 0.91;
    np3.add(hpturbine);

//	systemModule::ModulePtr recombinebleed2(new Mixer("recombinebleed2", np3.nodes_[17].get(), np3.nodes_[18].get(), np3.nodes_[8].get()));
//	Mixer *recombinebleed2C = (Mixer *) recombinebleed2.get();
//	recombinebleed2C->verbosity_ = 0;
//	recombinebleed2C->chokingFeedback_ = inlet;
//	recombinebleed2C->eff_ = 0.99;
//	np3.add(recombinebleed2);


	systemModule::ModulePtr link3(new EffIsenDuct("linkhpt2lpt", np3.nodes_[16].get(), np3.nodes_[17].get()));
	EffIsenDuct *link3C = (EffIsenDuct *) link3.get();
	link3C->etap_ = 0.995;
	link3C->etaT_ = 1.0;
	np3.add(link3);


    systemModule::ModulePtr lpturbine(new TurbineBleed("lpturbine", np3.nodes_[17].get(),  np3.nodes_[18].get(), np3.nodes_[8].get(), np3.mech_[0].get()));
    TurbineBleed *lpturbineC = (TurbineBleed *) lpturbine.get();
    lpturbineC->efft_ = 0.92;
	np3.add(lpturbine);

	systemModule::ModulePtr link4(new EffIsenDuct("linklpt2nozzle", np3.nodes_[18].get(), np3.nodes_[19].get()));
	EffIsenDuct *link4C = (EffIsenDuct *) link4.get();
	link4C->etap_ = 0.99;
	link4C->etaT_ = 1.0;
	np3.add(link4);

	systemModule::ModulePtr corenozzle(new AdaptedNozzle<EffIsenDuct>("corenozzle", np3.nodes_[19].get(),  np3.nodes_[20].get()));
	AdaptedNozzle<EffIsenDuct> *corenozzleC = (AdaptedNozzle<EffIsenDuct> *) corenozzle.get();
	corenozzleC->etap_ = 1.0;
	corenozzleC->chokingMach_ = 0.99;
	corenozzleC->freeStream_ =  &inlet->N1();
	np3.add(corenozzle);

	systemModule::ModulePtr link5(new EffIsenDuct("linkbypassduct", np3.nodes_[2].get(), np3.nodes_[21].get()));
	EffIsenDuct *link5C = (EffIsenDuct *) link5.get();
	link5C->etap_ = 0.985;
	link5C->etaT_ = 1.0;
	np3.add(link5);

	systemModule::ModulePtr bypassnozzle(new AdaptedNozzle<EffIsenDuct>("bypassnozzle", np3.nodes_[21].get(),  np3.nodes_[22].get()));
	AdaptedNozzle<EffIsenDuct> *bypassnozzleC = (AdaptedNozzle<EffIsenDuct> *) bypassnozzle.get();
	bypassnozzleC->etap_ = 1.0;
	bypassnozzleC->chokingMach_ = 0.99;
	bypassnozzleC->freeStream_ =  &inlet->N1();
	np3.add(bypassnozzle);

//    systemModule::ModulePtr recombine(new mixerModule("recombine", np3.nodes_[14].get(), &np3.N2(), np3.nodes_[16].get()));
//    mixerModule *recombineC = (mixerModule *) recombine.get();
//    recombineC->verbosity_ = 0;
//    recombineC->chokingFeedback_ = inlet;
//    recombineC->eff_ = 0.99;
//    np3.add(recombine);


// DONE UNDOING


  	np3.calculate();
	std::cout << std::endl <<"----------------------------"  << std::endl << "system.calculate() complete" << std::endl <<"----------------------------" << std::endl;
    returned.sys = new systemModule(np3);
//    returned.thrust = np3.thrust();

//    fmt::print(freeStream.report());
    double Tin = -freeStream.A()*freeStream.rho()*freeStream.getU()*freeStream.getU();

    double Toutb = bypassnozzle->N2().A()*bypassnozzle->N2().rho()*bypassnozzle->N2().getU()*bypassnozzle->N2().getU() - (bypassnozzle->N2().A()*(bypassnozzle->N2().getPress() - exit.getPress()));
    double Toutc = corenozzle->N2().A()*corenozzle->N2().rho()*corenozzle->N2().getU()*corenozzle->N2().getU() - (corenozzle->N2().A()*(corenozzle->N2().getPress() - exit.getPress()));


    returned.thrust = Toutb + Toutc + Tin;

    return {true, returned};


}

