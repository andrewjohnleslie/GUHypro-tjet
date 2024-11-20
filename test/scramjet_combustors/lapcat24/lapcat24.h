//
// Created by robert on 13/11/17.
//

#ifndef HYPRO_BURROWS_H
#define HYPRO_BURROWS_H

#endif // HYPRO_BURROWS_H

#include "core/Node.h"
#include "core/systemModule.h"

#include "InjectionPhi.h"
#include "ParMixer.h"
#include "balanceMach.h"
#include "combustion/EffComb.h"
#include "combustion/EquilCombustion.h"
#include "fmt/format.h"
#include "injectors/injection.h"
#include "simpleSplitter.h"
#include "solvers/EffIsenDuct.h"
#include "solvers/Friction.h"
#include "solvers/isentropicDuct.h"
// #include "Graphics/SystemGraph.h"

hypro::systemModule *burrows(
	const vector<double> A, const int combustionMode, const double inlet_p, const double inlet_T,
	const double inlet_U, const Cantera::compositionMap inletComposition,
	const double phi, const double combeff, const std::string chemistryModel)
{

	bool par_mixer = false;
	bool injector = true;
	bool premixed = false;

	// Create inlet node
	hypro::Node &inlet = *new hypro::Node(
		"/home/robert/repos/HyPro/test/"
		"scramjet_combustors/Burrows/gri30_highT.xml",
		chemistryModel);

	// Set inlet node thermofluiddynamic state
	inlet.Name_ = "Inlet";
	inlet.setTPY(inlet_T, inlet_p, inletComposition);
	inlet.setU(inlet_U);

	vector<double> empty_thermos(inlet.X().size());

	Cantera::compositionMap injectionComposition;
	injectionComposition.emplace("H2", 1);
	//   injection.setTPY(248.56, 92923.3596, injectionComposition);
	// Create exit node
	hypro::Node &exit = *new hypro::Node(inlet);
	exit.Name_ = "Exit";

	// Create system
	hypro::systemModule system("system", &inlet, &exit);

	// Create other nodes
	int numNodes = A.size();
	system.initNodes(numNodes - 2, inlet);

	// Assign node areas
	inlet.A(A[0]);
	for (int i = 0; i < numNodes - 2; i++)
	{
		system.nodes_[i]->A(A[i + 1]);
	}
	exit.A(A[numNodes - 1]);

	system.verbosity_ = 1;

	if (injector == true)
	{
		// ################################Injector
		hypro::systemModule::ModulePtr injector(
			new hypro::InjectionPhi<hypro::balanceMach>("injector", &system.N1(),
														system.nodes_[0].get()));
		injector->N2().Name_ = "Injection";
		hypro::InjectionPhi<hypro::balanceMach> *injectorC =
			(hypro::InjectionPhi<hypro::balanceMach> *)injector.get();
		injectorC->PhiMax_ = phi; // TODO Add throttle for injector here
		injectorC->Phi_ = phi;
		// secondaryC->mfrMax_ = *throttle;
		vector<double> Xf = empty_thermos;
		Xf[0] = 1;
		injectorC->Xf_ = Xf;
		injectorC->Tf_ = 248;
		injectorC->verbosity_ = 0;
		int rcoeff2;
		if (inlet.gasmodel_.get()->id() == "gri30")
		{
			rcoeff2 = 3;
		}
		else
		{
			rcoeff2 = 1;
		}

		injectorC->Rcoeff_[0] = 2;
		injectorC->Rcoeff_[rcoeff2] = 1;
		injector->chokingFeedback_ = injector;
		system.add(injector);

		hypro::systemModule::ModulePtr areachange2(new hypro::isentropicDuct(
			"initial change", system.nodes_[0].get(), system.nodes_[1].get()));
		areachange2->N2().Name_ = "duct end";
		system.add(areachange2);

		hypro::systemModule::ModulePtr isolator(new hypro::Friction(
			"Isolator", system.nodes_[1].get(), system.nodes_[2].get()));
		isolator->N2().Name_ = "Isolator end";
		hypro::Friction *isolatorC = (hypro::Friction *)isolator.get();
		isolatorC->L_ = 0.23;
		isolatorC->D_ = 2 * 0.093;										 // 4 * 0.026 * 1 / (2 * 1);
		isolatorC->Cf_ = 0.00070508377998784754 + 0.0032510047534616247; //isolatorFriction;
		isolatorC->heat_flux = true;
		isolatorC->qdot = 760959 - 32841;
		system.add(isolator);
	}

	if (par_mixer == true)
	{
		system.nodes_[0]->setTPY(248.56, 92923.3596, injectionComposition);
		system.nodes_[0]->setU(1284.9);
		hypro::systemModule::ModulePtr injector(new hypro::ParMixer("Injector", &system.N1(), system.nodes_[1].get(), system.nodes_[0].get()));
		injector->N2().Name_ = "Injector Out";
		system.add(injector);

		hypro::systemModule::ModulePtr areachange2(new hypro::isentropicDuct(
			"initial change", system.nodes_[1].get(), system.nodes_[2].get()));
		areachange2->N2().Name_ = "duct end";
		system.add(areachange2);

		hypro::systemModule::ModulePtr isolator(new hypro::Friction(
			"Isolator", system.nodes_[2].get(), system.nodes_[3].get()));
		isolator->N2().Name_ = "Isolator end";
		hypro::Friction *isolatorC = (hypro::Friction *)isolator.get();
		isolatorC->L_ = 0.23;
		isolatorC->D_ = 2 * 0.093;										 // 4 * 0.026 * 1 / (2 * 1);
		isolatorC->Cf_ = 0.00070508377998784754 + 0.0032510047534616247; //isolatorFriction;
		isolatorC->heat_flux = true;
		isolatorC->qdot = 760959 - 32841;
		system.add(isolator);
	}

	// ############################## Isolator
	//   hypro::systemModule::ModulePtr isolator(new hypro::Friction(
	//       "Isolator", system.nodes_[0].get(), system.nodes_[1].get()));

	if (combustionMode == 0)
	{
		// ====================== EffComb ======================
		hypro::systemModule::ModulePtr combustor(
			new hypro::EffComb<hypro::Friction>("combustor",
												system.nodes_[0].get(),
												system.nodes_[1].get()));
		combustor->N2().Name_ = "Chamber End";
		hypro::EffComb<hypro::Friction> *combustorC =
			(hypro::EffComb<hypro::Friction> *)combustor.get();
		if (inlet.gasmodel_.get()->id() == "gri30")
		{
			combustorC->Rcoeff_[0] = 2;
			combustorC->Rcoeff_[3] = 1;
			combustorC->Pcoeff_[5] = 2;
		}
		else
		{
			combustorC->Rcoeff_[0] = 2;
			combustorC->Rcoeff_[1] = 1;
			combustorC->Pcoeff_[3] = 2;
		}
		combustorC->eta_ = combeff;
		combustorC->verbosity_ = 0;
		// combustorC->chokingFeedback_ = injector;
		combustorC->Cf_ = 0.001264888478204313 + 0.0029981004934882564; //combustorFriction;
		combustorC->D_ = 2 * 0.093;										// 4*0.026*1/(2*1);
		combustorC->L_ = 0.126;
		combustorC->heat_flux = true;
		system.add(combustor);

		hypro::systemModule::ModulePtr areachange(new hypro::isentropicDuct(
			"areachange", system.nodes_[1].get(), &system.N2()));
		system.add(areachange);
	}
	else if (combustionMode == 1)
	{
		// ====================== EquilComb ======================
		hypro::systemModule::ModulePtr combustor(
			new hypro::EquilCombustion<hypro::Friction>(
				"combustor", system.nodes_[0].get(),
				system.nodes_[1].get()));
		combustor->N2().Name_ = "Chamber End";
		hypro::EquilCombustion<hypro::Friction> *combustorC =
			(hypro::EquilCombustion<hypro::Friction> *)combustor.get();
		combustorC->verbosity_ = 0;
		// combustorC->chokingFeedback_ = injector;
		combustorC->Cf_ = 0.001264888478204313 + 0.0029981004934882564; //combustorFriction;
		combustorC->D_ = 2 * 0.093;
		combustorC->L_ = 0.126;
		combustorC->heat_flux = true;
		system.add(combustor);

		hypro::systemModule::ModulePtr areachange(new hypro::isentropicDuct(
			"areachange", system.nodes_[1].get(), &system.N2()));
		system.add(areachange);
	}
	else if (combustionMode == 2)
	{
		// ====================== Split EffComb ======================
		std::vector<int> split_indexes;
		std::vector<int> isolator_indexes;
		std::vector<int> combustor_indexes;
		std::vector<int> mixer_indexes;
		std::vector<int> area_indexes;
		if (par_mixer == true)
		{
			split_indexes = {3, 4, 5};
			isolator_indexes = {4, 6};  // With parMixer
			combustor_indexes = {5, 7}; // With parMixer
			mixer_indexes = {6, 8, 7};  // With parMixer
			area_indexes = {8};			// With parMixer
		}
		else
		{
			split_indexes = {2, 3, 4};  // Without parMixer
			isolator_indexes = {3, 5};  // Without parMixer
			combustor_indexes = {4, 6}; // Without parMixer
			mixer_indexes = {5, 7, 6};  // Without parMixer
			area_indexes = {7};			// Without parMixer
		}

		std::cout << split_indexes[0] << " " << split_indexes[1] << " " << split_indexes[2] << std::endl;
		hypro::systemModule::ModulePtr splitter(new hypro::simpleSplitter(
			"splitter", system.nodes_[split_indexes[0]].get(),
			system.nodes_[split_indexes[1]].get(), system.nodes_[split_indexes[2]].get()));
		splitter->N2().Name_ = "Split 1";
		//        splitter->N3().Name_ = "Split 2";
		system.add(splitter);

		std::cout << isolator_indexes[0] << " " << isolator_indexes[1] << std::endl;
		hypro::systemModule::ModulePtr split_isolator(new hypro::Friction(
			"split_Isolator", system.nodes_[isolator_indexes[0]].get(), system.nodes_[isolator_indexes[1]].get()));
		split_isolator->N2().Name_ = "Split Isolator end";
		hypro::Friction *split_isolatorC = (hypro::Friction *)split_isolator.get();
		split_isolatorC->L_ = 0.126;
		split_isolatorC->D_ = 2 * 0.093 / 2;
		split_isolatorC->Cf_ = (0.001264888478204313 + 0.0029981004934882564) / 2;
		system.add(split_isolator);

		std::cout << combustor_indexes[0] << " " << combustor_indexes[1] << std::endl;
		hypro::systemModule::ModulePtr split_combustor(
			new hypro::EffComb<hypro::Friction>("split_combustor",
												system.nodes_[combustor_indexes[0]].get(),
												system.nodes_[combustor_indexes[1]].get()));
		split_combustor->N2().Name_ = "Split Chamber End";
		hypro::EffComb<hypro::Friction> *split_combustorC =
			(hypro::EffComb<hypro::Friction> *)split_combustor.get();
		if (inlet.gasmodel_.get()->id() == "gri30")
		{
			split_combustorC->Rcoeff_[0] = 2;
			split_combustorC->Rcoeff_[3] = 1;
			split_combustorC->Pcoeff_[5] = 2;
		}
		else
		{
			split_combustorC->Rcoeff_[0] = 2;
			split_combustorC->Rcoeff_[1] = 1;
			split_combustorC->Pcoeff_[3] = 2;
		}
		split_combustorC->eta_ = combeff;
		// split_combustorC->chokingFeedback_ = injector;
		split_combustorC->Cf_ = (0.001264888478204313 + 0.0029981004934882564) / 2;
		split_combustorC->D_ = 2 * 0.093 / 2;
		split_combustorC->L_ = 0.126;
		split_combustorC->heat_flux = true;
		split_combustorC->qdot = 760959 - 32841;
		system.add(split_combustor);

		std::cout << mixer_indexes[0] << " " << mixer_indexes[1] << " " << mixer_indexes[2] << std::endl;
		hypro::systemModule::ModulePtr mixer(
			new hypro::ParMixer("mixer", system.nodes_[mixer_indexes[0]].get(),
								system.nodes_[mixer_indexes[1]].get(), system.nodes_[mixer_indexes[2]].get()));
		system.add(mixer);

		std::cout << area_indexes[0] << std::endl;
		hypro::systemModule::ModulePtr areachange(new hypro::isentropicDuct(
			"areachange", system.nodes_[area_indexes[0]].get(), &system.N2()));
		system.add(areachange);
	}
	else if (combustionMode == 3)
	{
		// !!!!!!!!!!!!!!!!!!!!!!!! Split EquilComb !!!!!!!!!!!!!!!!!!!!!!!!!!
		std::vector<int> split_indexes;
		std::vector<int> isolator_indexes;
		std::vector<int> combustor_indexes;
		std::vector<int> mixer_indexes;
		std::vector<int> area_indexes;
		if (par_mixer == true)
		{
			split_indexes = {2, 3, 4};
			isolator_indexes = {3, 5};  // With parMixer
			combustor_indexes = {4, 6}; // With parMixer
			mixer_indexes = {5, 7, 6};  // With parMixer
			area_indexes = {7};			// With parMixer
		}
		else
		{
			split_indexes = {0, 1, 2};  // Without parMixer
			isolator_indexes = {1, 3};  // Without parMixer
			combustor_indexes = {2, 4}; // Without parMixer
			mixer_indexes = {3, 5, 4};  // Without parMixer
			area_indexes = {5};			// Without parMixer
		}

		hypro::systemModule::ModulePtr splitter(new hypro::simpleSplitter(
			"splitter", system.nodes_[split_indexes[0]].get(),
			system.nodes_[split_indexes[1]].get(), system.nodes_[split_indexes[2]].get()));
		splitter->N2().Name_ = "Split 1";
		//        splitter->N3().Name_ = "Split 2";
		system.add(splitter);

		hypro::systemModule::ModulePtr split_isolator(new hypro::Friction(
			"split_Isolator", system.nodes_[isolator_indexes[0]].get(), system.nodes_[isolator_indexes[1]].get()));
		split_isolator->N2().Name_ = "Split Isolator end";
		hypro::Friction *split_isolatorC = (hypro::Friction *)split_isolator.get();
		split_isolatorC->L_ = 0.126;
		split_isolatorC->D_ = 2 * 0.093 / 2;
		split_isolatorC->Cf_ = (0.001264888478204313 + 0.0029981004934882564) / 2;
		system.add(split_isolator);

		hypro::systemModule::ModulePtr split_combustor(
			new hypro::EquilCombustion<hypro::Friction>(
				"split_combustor", system.nodes_[combustor_indexes[0]].get(),
				system.nodes_[combustor_indexes[0]].get()));
		split_combustor->N2().Name_ = "Split Chamber End";
		hypro::EquilCombustion<hypro::Friction> *split_combustorC =
			(hypro::EquilCombustion<hypro::Friction> *)split_combustor.get();
		split_combustorC->verbosity_ = 0;
		// split_combustorC->chokingFeedback_ = injector;
		split_combustorC->Cf_ = (0.001264888478204313 + 0.0029981004934882564) / 2;
		split_combustorC->D_ = 2 * 0.093 / 2;
		split_combustorC->L_ = 0.126;
		split_combustorC->heat_flux = true;
		split_combustorC->qdot = 760959 - 32841;
		system.add(split_combustor);

		hypro::systemModule::ModulePtr mixer(
			new hypro::ParMixer("mixer", system.nodes_[mixer_indexes[0]].get(),
								system.nodes_[mixer_indexes[1]].get(), system.nodes_[mixer_indexes[2]].get()));
		system.add(mixer);

		hypro::systemModule::ModulePtr areachange(new hypro::isentropicDuct(
			"areachange", system.nodes_[area_indexes[0]].get(), &system.N2()));
		system.add(areachange);
	}

	system.calculate();
	hypro::Collection::OutType out = system.out();
	vector<double> propellant = system.propellantMfr();

	fmt::print("{}, {}, {}, {}, {}, {}, {}, ", 
		inlet.gasmodel_.get()->id(), system.N2().rho(), system.N2().getU(), 
		system.N2().getTemp(), system.N2().getPress(), system.N2().mfr(), 
		system.N2().T0());


	// std::cout << inlet.gasmodel_.get()->id() << ", " << system.N2().rho() << ", "
	// 		  << system.N2().getU() << ", " << system.N2().getTemp() << ", "
	// 		  << system.N2().getPress() << ", " << system.N2().mfr() << ", "
	// 		  << system.N2().T0() << ", ";
	Cantera::compositionMap outlet = system.N2().getYMap();
	// std::cout << outlet["H2"] << ", " << outlet["O2"] << ", " << outlet["N2"]
	// 		  << ", " << outlet["H2O"] << ", " << outlet["OH"] << ", "
	// 		  << outlet["O"] << ", " << outlet["H"] << ", " << outlet["NO"]
	// 		  << ", " << outlet["NO2"] << std::endl;
	fmt::print("{}, {}, {}, {}, {}\n", outlet["H2"], outlet["O2"], outlet["N2"],
		outlet["H2O"], outlet["OH"], outlet["O"], outlet["H"], outlet["NO"],
		outlet["NO2"]);		  
	// system.printOut(cout);
	std::cout << system.nodes_[1].get()->T0() << std::endl;
	// char* argv[0];
	// const int argc = sizeof(argv)/sizeof(char*);
	// system.show(0,0);

	// Save as JSON
	std::stringstream ss;
	system.saveJson(ss);

	// create and open a character archive for output
	std::ofstream ofs("burrows_injector");
	ofs << ss.str();
	ofs.close();

	return new hypro::systemModule(system);
}