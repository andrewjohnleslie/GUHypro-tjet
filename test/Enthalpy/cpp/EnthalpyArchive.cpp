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





//using namespace hypro;
//
//std::pair<bool, outputData>
//Execution(const int &verbosity) {
//
//   outputData returned;
//   int ondesign = 1.0;
//
//   double p, T, Mach;
//
////   T = 288.16;
//   T = 298.1;
//   p = 101325;
//   Mach = 0.01;
//
//   Node &freeStream =*new Node("/home/andrew/GUHypro-tjet/test/TurboJet/gri30_highT_modified_jetA.xml", "gri30");
//   freeStream.Name_ = "Freestream";
//
//   Node &exit = *new Node(freeStream);
//   exit.Name_ = "Nozzle exit to FS";
//
//   Cantera::compositionMap initialComposition;
//   initialComposition.emplace("O2", 0.209);
//   initialComposition.emplace("N2", 1 - 0.209);
//
//   freeStream.setTPX(T, p, initialComposition);
//   vector<double> empty_thermos(freeStream.X().size());
//   freeStream.setU(Mach * freeStream.geta());
//
//   systemModule system("system", &freeStream, &exit);
//   system.verbosity_ = verbosity;
//
//   double A = 5.725;
//
//   system.initNodes(4, freeStream);
//   freeStream.A(A);
//   system.nodes_[0]->A(A);
//   system.nodes_[1]->A(A);
//   system.nodes_[2]->A(A);
//   system.nodes_[3]->A(A);
//
//
//   exit.A(A);
//
//   system.initShaft(1);
//
//
////
////
////   systemModule::ModulePtr inlet(new EffIsenDuct("compressor", &system.N1(), system.nodes_[0].get()));
////   inlet->N2().Name_ = "Compressor Inlet";
////   system.add(inlet);
////
////   systemModule::ModulePtr compressor(new Compressor("compressor", system.nodes_[0].get(), system.nodes_[1].get(), system.mech_[0].get()));
////   compressor->N2().Name_ = "Compressor Exit";
////   Compressor *compressorC = (Compressor *) compressor.get();
////   compressorC->effc_ = 0.9;
////   compressorC->pi_c_ = 5.5;
////   compressorC->checkdesignpoint_ = ondesign;
////   system.add(compressor);
////
////   systemModule::ModulePtr injection(new InjectionPhi<balanceMach>("injection", system.nodes_[1].get(),  system.nodes_[2].get()));
////   injection->N2().Name_ = "Combustion Inlet";
////   InjectionPhi<balanceMach> *injectionC = (InjectionPhi<balanceMach> *) injection.get();
////   injectionC->Phi_ = 0.21775;
////   injectionC->PhiMax_ = 1.0;
////   injectionC->chokingFeedback_ = injection;
////   injectionC->chokingMach_ = 0.98;
////   injectionC->Tf_ = 495;
////   vector<double> Xf = empty_thermos;
////   Xf[53] = 1;
////   injectionC->Rcoeff_[3] = 71;
////   injectionC->Rcoeff_[53] = 4;
////   injectionC->Xf_ = Xf;
////   system.add(injection);
////
////
////
////
////	systemModule::ModulePtr combustion(new EffComb<balanceMach>("combustor", system.nodes_[2].get(),  system.nodes_[3].get()));
////	combustion->N2().Name_ = "Turbine Inlet";
////	EffComb<balanceMach> *combustionC = (EffComb<balanceMach> *) combustion.get();
////	combustionC->Rcoeff_[53] = 4;
////	combustionC->Rcoeff_[3]  = 71;
////	combustionC->Pcoeff_[15] = 48;
////	combustionC->Pcoeff_[5]  = 46;
////	combustionC->eta_ = 1.0;
////	combustionC->chokingFeedback_ = injection;
////	system.add(combustion);
////
////
////   systemModule::ModulePtr turbine(new Turbine("turbine", system.nodes_[3].get(), &system.N2(),  system.mech_[0].get()));
////   Turbine *turbineC = (Turbine *) turbine.get();
////   turbineC->efft_ = 0.85;
////   turbineC->checkdesignpoint_ = ondesign;
////   system.add(turbine);
////
////   system.calculate();
////   system.printOut(std::cout);
////
////   std::cout << freeStream.Name_ << "\n";
////   std::cout << freeStream.report();
////
////	for (std::size_t i = 0; i<(system.nodes_.size()); i++) {
////		   std::cout << system.nodes_[i]->Name_ << "\n";
////		   std::cout << system.nodes_[i]->report();
////	}
////
////
////   std::cout << exit.Name_ << "\n";
////   std::cout << exit.report();
//
//
////   Node &Dummy = *new Node(freeStream);
////   Dummy.Name_ = "Dummy";
////   Dummy.A(A);
////   Dummy.setTP(298.1115,p);
////   Dummy.setU(0);
//
//	Node &Test = *new Node(freeStream);
//	Test.Name_ = "Test";
//	Test.A(1);
//	Test.setTP(350,p);
//	Test.setU(10);
//
//	Node &TestOut = *new Node(freeStream);
//	TestOut.Name_ = "TestOut";
//	TestOut.A(1);
//
//    systemModule testsystem("testsystem", &Test, &TestOut);
//    testsystem.initShaft(1);
//
//
//	systemModule::ModulePtr turbine(new Turbine("turbine", &testsystem.N1(), &testsystem.N2(),  testsystem.mech_[0].get()));
//	Turbine *turbineC = (Turbine *) turbine.get();
//	turbineC->efft_ = 1.0;
//	testsystem.mech_[0]->setPower(100000);
//	turbineC->checkdesignpoint_ = ondesign;
//	testsystem.add(turbine);
//
//	std::cout << "\n\n" << "Beginning Test Node HyPro simulation " << "\n\n" << "************************************" << "\n\n";
//
//	testsystem.calculate();
//	testsystem.printOut(std::cout);
//
//
//	std::cout << "\n\nGas Model Enthalpy Mass: \t\t"<< Test.gasmodel_->enthalpy_mass() << "\n";
//
//
//
//
//	double Ttest, cptest;
//	Ttest = Test.getT0();
//	cptest = Test.Cp();
//	double Enth = cptest*Ttest;
//	std::cout << "Test Enthalpy (cp*T0) is: \t\t" << Enth << "\n";
//	std::cout << "Test Node MFR is: \t\t\t" << Test.mfr() << "\n";
//	std::cout << "Test node Density is: \t\t\t" << Test.rho() << "\n";
//	std::cout << "Test Enthalpy (.H()) is: \t\t" << Test.H() << " (Total) J/kg \n";
//	std::cout << "Test Enthalpy (.h()) is: \t\t" << Test.h() << " (Static) J/kg \n\n";
//
//	std::cout << "TestOut Node MFR is: \t\t\t" << TestOut.mfr() << "\n";
//	std::cout << "TestOut node Density is: \t\t" << TestOut.rho() << "\n";
//	std::cout << "TestOut Enthalpy (.H()) is: \t\t" << TestOut.H() << " (Total) J/kg \n";
//	std::cout << "TestOut Enthalpy (.h()) is: \t\t" << TestOut.h() << " (Static) J/kg \n\n";
//
//	std::cout << "\nDeltaH between input-output is: \t" << Test.H() - TestOut.H() << " (Total) J/kg \n\n";
//
//	//   std::cout << Test.report();
//
//
//	//   Dummy.~Node();
//	Test.~Node();
//	return {true, returned};
//
//
//}




using namespace Cantera;

void thermo_demo(double T, const std::string& file, const std::string& phase)
{


    shared_ptr<ThermoPhase> gas(newPhase(file,phase));
//    std::cout << typeid(*(gas)).name() << '\n';
//    gas->setState_TPX(T, OneAtm, "O2:0.209, H2:0.791");
//    gas->setState_TPX(T, OneAtm, "O2:0.161115, N2:0.779508, CO2:0.0303202, H2O:0.0290568");
//    gas->setState_TPX(T, OneAtm, "O2:0.151123, N2:0.77711, CO2:0.0366471, H2O:0.0351201");

    double PyCyclePressure = 101324.703484;
    double PyCyclepBurnerPressure = 1326848.29868;

    double HyPropBurnerPressure =  1316590;


//	Air Values
//    gas->setState_TPX(T, PyCyclePressure, "O2:0.209476, N2:0.78084, CO2:0.000319, Ar:0.009365");
//    gas->setState_TPX(T, PyCyclePressure,  "Ar:0.009365005078738, CO:0.000000002896520, CO2:0.000000002896520, N:0.000000002896520, NO:0.000000002896520, NO2:0.000000002896520, N2:0.780840504144310 ,O:0.000000002896520, O2:0.209476127566503");

//	Air-Fuel Values
//    gas->setState_TPX(T, HyPropBurnerPressure,  "Ar:0.00920267, H2O:0.0346681, CO2:0.0364888, N2:0.767305, O2:0.152336"); // HyPro Comp
    gas->setState_TPX(T, PyCyclePressure,  "Ar:0.009202688573271, H2O:0.034668606853779, CO2:0.036491242646628, NO:0.000369915115618, NO2:0.000011587004405, N2:0.767114150507551, OH:0.000003765776432, O2:0.152134181406649");  //pycyle comp
//    gas->setState_TPX(T, PyCyclePressure,  "C12H23:1.0");



//	gas->setState_TPX(T, OneAtm, "O2:1, H2:1, CO2:1");
//	gas->setState_TPX(T, OneAtm, "O2:0.5, H2:0.5");
//  gas->setState_TPX(T, OneAtm, "N2:1");


//    std::cout << "Beginning Cantera Demonstration" << "\n\n" << "*******************************" << "\n" << std::endl;


//    // temperature, pressure, and density
//    std::cout << "Gas Temperature is:\t" 	<< gas->temperature() 	<< "\n";
//    std::cout << "Gas Pressure is:\t" 		<< gas->pressure() 		<< "\n";
//    std::cout << "Gas Gamma is: \t" 	  	<< gas->cp_mass()/gas->cv_mass()		<< "\n\n";
//
//    // molar thermodynamic properties
//    std::cout << "Gas Enthalpy/mol is:\t" 	<< gas->enthalpy_mole() << "\n";
//    std::cout << "Gas Cp/mol is: \t\t" 		<< gas->cp_mole() 		<< "\n";
//
//    // specific (per unit mass) thermodynamic properties
//    std::cout << "Gas Enthalpy/kg is:\t"  	<< gas->enthalpy_mass() << "\n";
//    std::cout << "Gas Cp/kg is:\t\t"   		<< gas->cp_mass() 		<< "\n";



    std::cout << gas->report();
    // chemical potentials of the species
    int numSpecies = gas->nSpecies();
    double AMU = gas->meanMolecularWeight();
    vector_fp mu(numSpecies);
    vector_fp mucp(numSpecies);
    gas->getPartialMolarEnthalpies(&mu[0]);
    gas->getPartialMolarCp(&mucp[0]);
    vector_fp Xval(numSpecies);
    vector_fp Wval(numSpecies);
    gas->getMoleFractions(&Xval[0]);
    gas->getMolecularWeights(&Wval[0]);
    int n;
    double Hval {0};

    std::cout << "\n\n";
	fmt::print(std::cout, "{:<12} {:<12}\n", "Species", "Enthalpy");

    for (n = 0; n < numSpecies; n++) {

    	if (Xval[n] != 0)

    	{
//    		std::cout << gas->speciesName(n) << "MolarEnthalpy:\t\t " << ((mu[n])) << std::endl;
//    		std::cout << gas->speciesName(n) << "MolarEnthalpy*MolFraction:\t " << ((mu[n])) << std::endl;
    		string gasstring = (gas->speciesName(n)) + ':';
    		fmt::print(std::cout, "{:<12} {:<12.9f}\n", (gasstring), ((mu[n])));

    		Hval += (Xval[n]*mu[n]);
    	}

    }
    std::cout << "\nEnthalpy/mol is: "<<(Hval)<<'\n';


//    hypro::Node &freeStream =*new Node(file, "gri30");
//    compositionMap initialComposition;
//    initialComposition.emplace("O2", 0.161115);
//    initialComposition.emplace("N2", 0.779508);
//    initialComposition.emplace("CO2", 0.0303202);
//    initialComposition.emplace("H2O", 0.0290568);
//    freeStream.setTPX(350.0, 101325, initialComposition);
////    vector<double> PlaceHolder = freeStream.mfrX();
////    int sizePlace = PlaceHolder.size();
//    double PlaceHolder = freeStream.Nfr();
//    std::cout << "Nfr: \t\t" << PlaceHolder << "\n";


//    std::cout <<  gas->enthalpy_mole() << std::endl;


}


//void simple_demo2()
//{
//    // Create a new phase
////    std::unique_ptr<ThermoPhase> gas(newPhase("gri30.cti", "gri30_mix"));
//    std::unique_ptr<ThermoPhase> gas(newPhase("/home/andrew/GUHypro-tjet/test/Enthalpy/gri30_highT_modified_jetA.xml", "gri30_mix"));
//
//    // List of phases participating in reactions (just one for homogeneous
//    // kinetics)
//    std::vector<ThermoPhase*> phases{gas.get()};
//
//    // Create the Kinetics object. Based on the phase definition used, this will
//    // be a GasKinetics object.
//    std::unique_ptr<Kinetics> kin(newKineticsMgr(gas->xml(), phases));
//
//    // Set an "interesting" mixture state where we will observe non-zero reacton
//    // rates.
//    gas->setState_TPX(500.0, 2.0*OneAtm, "C12H23:1.0, O2:1.0, N2:3.76");
//    gas->equilibrate("HP");
//    gas->setState_TP(gas->temperature() - 100, gas->pressure());
//
//    // Get the net reaction rates
//    vector_fp wdot(kin->nReactions());
//    kin->getNetRatesOfProgress(wdot.data());
//
//    writelog("Net reaction rates for reactions involving C12H23\n");
//    size_t kCO2 = gas->speciesIndex("O2");
//    for (size_t i = 0; i < kin->nReactions(); i++) {
//        if (kin->reactantStoichCoeff(kCO2, i)
//            || kin->productStoichCoeff(kCO2, i)) {
//            writelog("{:3d}  {:30s}  {: .8e}\n",
//                i, kin->reactionString(i), wdot[i]);
//        }
//    }
//    writelog("\n");
//
//    // Create a Transport object. Based on the transport model specified in the
//    // "gri30_mix" phase, this will be a MixGasTransport object.
//    std::unique_ptr<Transport> trans(newDefaultTransportMgr(gas.get()));
//    writelog("T        viscosity     thermal conductivity\n");
//    writelog("------   -----------   --------------------\n");
//    for (size_t n = 0; n < 5; n++) {
//        double T = 300 + 100 * n;
//        gas->setState_TP(T, gas->pressure());
//        writelog("{:.1f}    {:.4e}    {:.4e}\n",
//            T, trans->viscosity(), trans->thermalConductivity());
//    }
//}

using namespace hypro;


void enthalpy_expansion(){


		outputData returned;
		double p, T, Mach;


//
//		T = 843.9885;
//		p = 820727.9085;
//		Mach = 0.28703022;
//		double A = 0.0236;

//		T = 1337.03;
//		p =  798632.468;
//		Mach = 0.252147;
//		double A = 0.319345;

		Node &freeStream =*new Node("/home/andrew/GUHypro-tjet/test/TurboJet/gri30_highT_modified_jetA.yaml", "gri30");
		freeStream.Name_ = "Freestream";

		T = 288.15;
		p =  101324.703484;
		Mach = 0.00001;
		double A = 1000;

	    Cantera::compositionMap initialComposition;
	    initialComposition.emplace("O2", 0.209476);
	    initialComposition.emplace("N2", 0.78084);
	    initialComposition.emplace("Ar", 0.009365);
	    initialComposition.emplace("CO2", 0.000319);

	    freeStream.setTPX(T, p, initialComposition);
		vector<double> empty_thermos(freeStream.X().size());
		freeStream.setU(0.00034138);


//		//	Bypass Composition
//		Cantera::compositionMap initialComposition;
//		initialComposition.emplace("O2", 0.209);
//		initialComposition.emplace("N2", 1 - 0.209);


////			Core Composition
//		Cantera::compositionMap initialComposition;
//		initialComposition.emplace("O2", 0.164749005043437);
//		initialComposition.emplace("N2", 0.780380158148393);
//		initialComposition.emplace("CO2", 0.028019150710555);
//		initialComposition.emplace("H2O", 0.0268516860976152);


//		freeStream.setTPX(T, p, initialComposition);
//		vector<double> empty_thermos(freeStream.X().size());
//		freeStream.setU(Mach * freeStream.geta());

		Node &exit = *new Node(freeStream);
		exit.Name_ = "Nozzle exit to FS";

		systemModule system("system", &freeStream, &exit);
		system.verbosity_ = 1;

//		exit.A(0.0284);


//		system.initNodes(1, freeStream);
//		system.nodes_[0]->A(A);

		double Ascale = 2.5;
		freeStream.A(A);
		exit.A(Ascale*A);

		system.initShaft(1);



//		systemModule::ModulePtr compressor(new Compressor("compressor", &system.N1(), &system.N2(), system.mech_[0].get()));
//		Compressor *compressorC = (Compressor *) compressor.get();
//		compressorC->effc_ = 1.00;
//		compressorC->pi_c_ = 2.6230;
//		system.add(compressor);




		systemModule::ModulePtr inlet(new EffIsenDuct("compressor", &system.N1() , &system.N2()));
		EffIsenDuct *inletC = (EffIsenDuct *) inlet.get();
		inlet->N1().Name_ = "Bleed Entrance";
		inlet->N2().Name_ = "Bleed Exit";
		system.add(inlet);





//		systemModule::ModulePtr turbine(new Turbine("turbine", &system.N1() , &system.N2(), system.mech_[0].get()));
//		Turbine *turbineC = (Turbine *) turbine.get();
//		turbine->N1().Name_ = "Bleed Entrance";
//		turbine->N2().Name_ = "Bleed Exit";
//		turbineC->efft_ = 1.0;
//		system.mech_[0]->setPower(210200*13.0873); // Bypass
////		system.mech_[0]->setPower(322500*118.389); // Core
////		system.mech_[0]->setPower(280000*118.389);
//		system.add(turbine);






////		double p45_ = 313105;
//		double p45_ = 306297; //Bypass Stream
////		double p45_ = 233925; //Core Stream
//		exit.isentropicP(freeStream, p45_);
//		double U2_ = (1/Ascale)*freeStream.mfr()/(A*exit.rho());
//		exit.setU(U2_);


		system.calculate();
		std::cout << "H is: " << system.N1().H() << "\n" ;
		std::cout << "h is: " << system.N1().h() << "\n" ;

  		system.printOut(std::cout);
  		std::cout << inletC->N1().report() << '\n';
//  		std::cout<<"test\n";



}


int main(int argc, char** argv) {


	enthalpy_expansion();
	return 0;

}


//int main()
//{
//    try {
//        simple_demo2();
//    } catch (std::exception& err) {
//        std::cout << err.what() << std::endl;
//    }
//}





//int main(int argc, char** argv) {
//
//	std::vector<double> T {1308.97};//{200, 300, 400, 500, 600, 700, 800, 900, 1000};
//
//	for (std::size_t Ti = 0; Ti < (T.size()); Ti++) {
//
//		std::cout << "\nT is: \t" << T[Ti] << "\n\n";
//
//		try {
//
////			thermo_demo(T[Ti], "h2o2.cti","ohmech");
//			thermo_demo(T[Ti], "/home/mark/git/GUHypro/test/Enthalpy/gri30_highT_modified_jetA.xml", "gri30_mix");
////			thermo_demo(T[Ti], "/usr/lib/python2.7/dist-packages/cantera/data/airNASA9.xml", "airNASA9");
//
//
//		}
//		catch (CanteraError& err) {
//
//			std::cout << err.what() << std::endl;
//			return 1;
//		}
//
//	}
//
//	std::cout << "\nExiting Main\n";
//	int verbosity = 1.0;
////	Execution(verbosity);
//
//
//
//    return 0;
//
//}


// NEW ADDED March 5th 2021


//void enthalpy_expansion(){
//
//
//		outputData returned;
//		double p, T, Mach;
//
//
//
//
////		Node &freeStream =*new Node("/home/mark/git/GUHypro/test/Enthalpy/gri30_highT_modified_jetA.xml", "gri30");
////		Node &freeStream =*new Node("/home/mark/git/GUHypro/test/Enthalpy/nasa.xml", "mark");
//		Node &freeStream =*new Node("/home/mark/git/GUHypro/test/Enthalpy/airNASA9.xml", "airNASA9");
//
//		freeStream.Name_ = "Freestream";
//
////		T = 288.15;
////		p =  101325.39296;
////		double A = 1.0000000;
////		Mach = 0.00001;
//
//
////		T = 656.364;
////		p =  1331072.97;
////		double A = 0.0933036424246;
////		Mach = 0.19995;
//
////		T = 1308.9854;
////		p =  1316548.48 ;
////		double A = 0.14152630920488;
////		Mach = 0.19981;
//
////		T = 1308.74605388889;
////		p =  1292825.00 ;
////		double A = 0.14152630920488;
////		Mach = 0.19981;
//
////		double T02 = 1316.957207;
////		double p02 =  1351194.92 ;
////		double Mach02 = 0.00000000;
//
//		T = 978.005356 ;
////		T = 1004.41809944444;
////		p =  307386.00;
////		T = 979.026491666667;
//		p =  508013.819;
//		double A = 0.2547292389276;
//		Mach = 0.40055;
//
//	    Cantera::compositionMap initialComposition;
////	    initialComposition.emplace("O2", 0.209476);
////	    initialComposition.emplace("N2", 0.78084);
////	    initialComposition.emplace("Ar", 0.009365);
////	    initialComposition.emplace("CO2", 0.000319);
//
//
//
//		// Turbine In Composition
//		initialComposition.emplace("O2", 0.152134831117542);
//		initialComposition.emplace("N2", 0.767116912136294);
//		initialComposition.emplace("Ar", 0.00920270079087897);
//		initialComposition.emplace("CO2", 0.0364914697290112);
//		initialComposition.emplace("NO2", 1.16187917087244e-05);
//		initialComposition.emplace("NO", 0.000370007095074531);
//		initialComposition.emplace("OH", 3.6293843787906e-06);
//		initialComposition.emplace("H2O", 0.0346687377373258);
//
//
////	    // Turbine Out Composition
////	    initialComposition.emplace("O2", 0.152315752009611);
////		initialComposition.emplace("N2", 0.767289899340093);
////		initialComposition.emplace("Ar", 0.00920266602059696);
////		initialComposition.emplace("CO2", 0.036491334525617);
////		initialComposition.emplace("NO2", 2.26350368790788e-06);
////		initialComposition.emplace("NO", 2.7588240861387e-05);
////		initialComposition.emplace("OH", 5.02663151113896e-08);
////		initialComposition.emplace("H2O", 0.0346704142289617);
//
//
//
//
//	    freeStream.setTPX(T, p, initialComposition);
//		vector<double> empty_thermos(freeStream.X().size());
//		freeStream.setU(Mach*freeStream.geta());
//
//
//
//
//
//		Node &exit = *new Node(freeStream);
//		exit.Name_ = "Nozzle exit to FS";
//
////	    exit.setTPX(T02, p02, initialComposition);
////		exit.setU(Mach02*freeStream.geta());
//
//
//		systemModule system("system", &freeStream, &exit);
//		system.verbosity_ = 1;
//		freeStream.A(A);
//		exit.A(0.25472923957276);
//		system.initShaft(1);
//
//
//
//	    systemModule::ModulePtr turbine(new Turbine("turbine", &system.N1(),  &system.N2(), system.mech_[0].get()));
//	    turbine->N2().Name_ = "Turbine Exit";
//	    Turbine *turbineC = (Turbine *) turbine.get();
//	    turbineC->efft_ = 0.839031;
//	    system.mech_[0]->setPower(2.568682939248800e+07);
//	    turbineC->checkdesignpoint_ = 1;
//		system.add(turbine);
//
//
//
//		std::cout << freeStream.report() ;
//		std::cout << "\n\n";
//		freeStream.printOut(std::cout,true);
////		exit.printOut(std::cout,true);
//
////		system.calculate();
//
////  		system.printOut(std::cout);
//
//
//		int numSpecies = freeStream.gasmodel_->nSpecies();
//		double AMU = freeStream.gasmodel_->meanMolecularWeight();
//		Cantera::vector_fp mu(numSpecies);
//		Cantera::vector_fp mus(numSpecies);
//		Cantera::vector_fp mucp(numSpecies);
//		freeStream.gasmodel_->getPartialMolarEnthalpies(&mu[0]);
//		freeStream.gasmodel_->getPartialMolarEntropies(&mus[0]);
////		freeStream.gasmodel_->getPartialMolarCp(&mucp[0]);
//		freeStream.gasmodel_->getEnthalpy_RT(&mucp[0]);
//		Cantera::vector_fp Xval(numSpecies);
//		Cantera::vector_fp Wval(numSpecies);
//		freeStream.gasmodel_->getMoleFractions(&Xval[0]);
//		freeStream.gasmodel_->getMolecularWeights(&Wval[0]);
//		double Hval, sval, cpval {0};
//
//		std::cout << "\n\n";
//		fmt::print(std::cout, "{:<20} {:<20} {:<20} {:<20} {:<20} {:<20} \n\n", "Species", "EnthalpyPMol", "H/RT",  "MassFraction",  "EnthalpyXMolFracPkmol", "EnthalpyXMolFracPkg");
//
//		for (int n = 0; n < numSpecies; n++) {
////
//			if (Xval[n] != 0)
//
//			{
//				string gasstring = (freeStream.gasmodel_->speciesName(n)) + ':';
//				fmt::print(std::cout, "{:<20} {:<20.4f} {:<20.10f} {:<20.12f} {:<20.4f} {:<20.4f}\n", (gasstring), ((mu[n])), (mucp[n]),  Xval[n], (Xval[n]*mu[n]), (Xval[n]*mu[n]/AMU));
////				fmt::print(std::cout, "{:<20} {:<20.9f} {:<20.9f} {:<20.9f} {:<20.9f} {:<20.9f} {:<20.9f} {:<20.9f}\n", (gasstring), ((mu[n])), (mucp[n]), ((mus[n])), Xval[n], (Xval[n]*mu[n]), (Xval[n]*mucp[n]), (Xval[n]*mus[n]));
//
////
//				Hval += (Xval[n]*mu[n]/AMU);
//				sval += (Xval[n]*mus[n]/AMU);
//				cpval += (Xval[n]*mucp[n]/AMU);
//
//			}
////
//		}
//
////		double RT = freeStream.RT();
////		double R = freeStream.gasmodel_->RT();
////
////		fmt::print(std::cout, "\n\n{:<20} {:<20.9f}", "RT is: ", RT);
////		fmt::print(std::cout, "\n{:<20} {:<20.9f}", "R is: ", R);
////		fmt::print(std::cout, "\n{:<20} {:<20.9f}", "Cantera Gas Constant Is: " , Cantera::GasConstant );
//
//		std::cout << "\nEnthalpy/kg is: "<<(Hval)<<'\n';
//
//}




