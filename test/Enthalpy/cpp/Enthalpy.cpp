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
#include "Mixer.h"
#include "thermo/thermoState.h"
#include "NeutralLink.h"


#include "cantera/thermo/ThermoPhase.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/Phase.h"
#include "cantera/kinetics.h"
#include "cantera/transport.h"
#include <memory>
#include "typedefs.h"
#include "offDesignData.h"






using namespace hypro;


void enthalpy_expansion(){


		outputData returned;
		double p, T, Mach;




//		Node &freeStream =*new Node("/home/andrew/GUHypro-tjet/test/Enthalpy/gri30_highT_modified_jetA.xml", "gri30");
//		Node &freeStream =*new Node("/home/mark/git/GUHypro/test/Enthalpy/nasa.xml", "mark");

//		Node &freeStream =*new Node("/home/mark/git/GUHypro/test/Enthalpy/airNASA9.xml", "airNASA9");
		Node &freeStream =*new Node("/home/andrew/GUHypro-tjet/test/Enthalpy/airNASA9.yaml", "airNASA9");		

		freeStream.Name_ = "Freestream";

		// Compressor Inlet
//		T = 268.775324;
//		p =  79432.07;
//		double A = 0.3298;
//		Mach = 0.59999;

		T = 367.84;
		p =  60047;
		double A = 1.0946;
		Mach = 0.549;


	    Cantera::compositionMap initialComposition;
//	    initialComposition.emplace("O2", 0.209476);
//	    initialComposition.emplace("N2", 0.78084);
//	    initialComposition.emplace("Ar", 0.009365);
//	    initialComposition.emplace("CO2", 0.000319);
	    initialComposition.emplace("O2", 0.209);
	    initialComposition.emplace("N2", 1 - 0.209);






	    freeStream.setTPX(T, p, initialComposition);
		vector<double> empty_thermos(freeStream.X().size());
		freeStream.setU(Mach*freeStream.geta());





		Node &exit = *new Node(freeStream);
		exit.Name_ = "Nozzle exit to FS";

//	    exit.setTPX(T02, p02, initialComposition);
//		exit.setU(Mach02*freeStream.geta());


		systemModule system("system", &freeStream, &exit);
		system.verbosity_ = 1;
		freeStream.A(A);
		exit.A(0.0933);
		system.initShaft(1);



	    systemModule::ModulePtr compressor(new Compressor("compressor", &system.N1(),  &system.N2(), system.mech_[0].get()));
	    compressor->N2().Name_ = "Compressor Exit";
	    Compressor *compressorC = (Compressor *) compressor.get();
	    system.mech_[0]->eff_mech_ = 1.0;
		compressorC->effc_ = 0.87830073;
		compressorC->pi_c_ = 13.500000000;
		compressorC->checkdesignpoint_ = 1;
//	        if (ondesign == 0){
//	    			compressorC->setT04(T04_);
//	    			compressorC->kH_ = kH_;
		system.mech_[0]->Cp2_ =  1228;
		system.mech_[0]->Cp_ =  1228;
		system.mech_[0]->Ht_ =  24000000;
		system.mech_[0]->mfr_ =  freeStream.mfr();
//		system.mech_[0]->kH_ = kH_;
//	    			compressorC->effc_ = 0.87830073;
//	    		}
		system.add(compressor);




		std::cout << freeStream.h();
		std::cout << "\n\n";
		freeStream.printOut(std::cout,true);
//		exit.printOut(std::cout,true);

//		system.calculate();
//		double H02 = 369231.641 ;
//		compressorC->PressRise(H02);
//  		system.printOut(std::cout);



}


int main(int argc, char** argv) {


	enthalpy_expansion();
	return 0;

}







