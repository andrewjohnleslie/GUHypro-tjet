//
// Created by robert on 13/11/17.
//

#ifndef HYPRO_BURROWS_H
#define HYPRO_BURROWS_H

#endif //HYPRO_BURROWS_H

#include "core/Node.h"
#include "core/systemModule.h"

#include "InjectionPhi.h"
#include "injectors/injection.h"
#include "balanceMach.h"
#include "combustion/EffComb.h"
#include "solvers/Friction.h"
#include "combustion/EquilCombustion.h"
#include "combustion/separatedFlowCombustor.h"
#include "ParMixer.h"
#include "simpleSplitter.h"
#include "solvers/isentropicDuct.h"
#include "solvers/EffIsenDuct.h"

hypro::systemModule* lorrain(const int numNodes, const vector<double> A, const double combustorFriction, const double combustorLength, const int combustorNode1, const int combustionMode, const double isolatorFriction, const double isolatorLength, const double inlet_p, const double inlet_T, const double inlet_U, const Cantera::compositionMap inletComposition, const double phi, const double combeff, const int variableDuct) {

   hypro::Node &inlet = *new hypro::Node("/home/robert/repos/Cantera_HyPro/HyPro/test/scramjet_combustors/Burrows/gri30_highT.xml", "gri30");
//    hypro::Node &inlet = *new hypro::Node("/home/robert/repos/Cantera_HyPro/HyPro/test/scramjet_combustors/Burrows/gri30_highT.xml", "edm");
//    hypro::Node &inlet = *new hypro::Node("/home/robert/repos/Cantera_HyPro/HyPro/test/scramjet_combustors/Burrows/gri30_highT.xml", "frc");
    inlet.Name_ = "Combustor_Inlet";
    inlet.setTPY(inlet_T, inlet_p, inletComposition);
    inlet.setU(inlet_U);
//    std::cout << inlet.report() << std::endl;

    vector<double> empty_thermos(inlet.X().size());

    hypro::Node &exit = *new hypro::Node(inlet);
    exit.Name_ = "Nozzle Exit";

    hypro::systemModule system("system", &inlet, &exit);
    system.initNodes(numNodes - 2,inlet);
    for (int i = 0; i < numNodes - 2; i++) {
        system.nodes_[i]->A(A[i + 1]);
    }
    inlet.A(A[1]);
    exit.A(A[numNodes - 1]);

    system.verbosity_ = 1;

//    // ################################ Injector
//    hypro::systemModule::ModulePtr injector(
//            new hypro::InjectionPhi<hypro::balanceMach>("injector", &system.N1(), system.nodes_[0].get()));
//    injector->N2().Name_ = "Injection";
//    hypro::InjectionPhi<hypro::balanceMach> *injectorC = (hypro::InjectionPhi<hypro::balanceMach> *) injector.get();
//    injectorC->PhiMax_ = phi; //TODO Add throttle for injector here
//    injectorC->Phi_ = phi;
//    //secondaryC->mfrMax_ = *throttle;
//    vector<double> Xf = empty_thermos;
//    Xf[0] = 1;
//    injectorC->Xf_ = Xf;
//    injectorC->Tf_ = 248;
//    injectorC->verbosity_ = 0;
//    int rcoeff2;
//    if (inlet.gasmodel_.get()->id() == "gri30"){
//        rcoeff2 = 3;
//    } else {
//        rcoeff2 = 1;
//    }
//
//    injectorC->Rcoeff_[0] = 2;
//    injectorC->Rcoeff_[rcoeff2] = 1;
//    injector->chokingFeedback_ = injector;
//    system.add(injector);

//    hypro::systemModule::ModulePtr injector(new hypro::injection(
//            "Injection",&system.N1(),system.nodes_[0].get()));
//    injector->N2().Name_ = "Injection end";
//    hypro::injection* injectorC = (hypro::injection*)injector.get();
//    vector<double> Xf = empty_thermos;
//    Xf[0] = 0.25;
//    injectorC->Xf_ = Xf;
//    system.add(injector);

    // ############################## Isolator
    hypro::systemModule::ModulePtr isolator(new hypro::Friction(
            "Isolator", &system.N1(), system.nodes_[0].get()));
    isolator->N2().Name_ = "Isolator end";
    hypro::Friction *isolatorC = (hypro::Friction *) isolator.get();
    isolatorC->L_ = isolatorLength;
    isolatorC->D_ = 4*0.026*1/(2*1); //4 * 0.026 * 1 / (2 * 1);
    isolatorC->Cf_ = isolatorFriction;
    system.add(isolator);
//    int nextNode = 1;

    if (combustionMode == 0){
        // ====================== EffComb ======================
        hypro::systemModule::ModulePtr combustor(new hypro::EffComb<hypro::Friction>("combustor",system.nodes_[0].get(),&system.N2()));
//        hypro::systemModule::ModulePtr combustor(new hypro::EffComb<hypro::Friction>("combustor",system.nodes_[0 + nextNode].get(),system.nodes_[1 + nextNode].get()));
//        hypro::systemModule::ModulePtr combustor(new hypro::EffComb<hypro::Friction>("combustor",system.nodes_[0 + nextNode].get(),system.nodes_[1 + nextNode].get()));
        combustor->N2().Name_ = "Chamber End";
        hypro::EffComb<hypro::Friction>* combustorC = (hypro::EffComb<hypro::Friction>*)combustor.get();
////        Full GRI30
        combustorC->Rcoeff_[0] = 2;
        combustorC->Rcoeff_[3] = 1;
        combustorC->Pcoeff_[5] = 2;
//        Reduced Species
//        combustorC->Rcoeff_[0] = 2;
//        combustorC->Rcoeff_[1] = 1;
//        combustorC->Pcoeff_[3] = 2;
        combustorC->eta_ = combeff;
        combustorC->verbosity_ = 0;
//        combustorC->chokingFeedback_ = injector;
        combustorC->Cf_ = combustorFriction;
        combustorC->D_ = 4*0.026*1/(2*1);//4*0.026*1/(2*1);
        combustorC->L_ = combustorLength;
        combustorC->heat_flux = true;
        system.add(combustor);

    } else if (combustionMode == 1){
//        // ====================== EquilComb ======================
//        hypro::systemModule::ModulePtr combustor(new hypro::EquilCombustion<hypro::Friction>("combustor",system.nodes_[0 + nextNode].get(),system.nodes_[1 + nextNode].get()));
//        combustor->N2().Name_ = "Chamber End";
//        hypro::EquilCombustion<hypro::Friction>* combustorC = (hypro::EquilCombustion<hypro::Friction>*)combustor.get();
//        combustorC->verbosity_ = 0;
//        combustorC->chokingFeedback_ = injector;
//        combustorC->Cf_ = combustorFriction;
//        combustorC->D_ = 2*0.089;
//        combustorC->L_ = combustorLength;
//        system.add(combustor);

    } else if (combustionMode == 2){
        // ====================== Split EffComb ======================
//        hypro::systemModule::ModulePtr splitter(new hypro::simpleSplitter("splitter",system.nodes_[0 + nextNode].get(),system.nodes_[1 + nextNode].get(),system.nodes_[2 + nextNode].get()));
//        splitter->N2().Name_ = "Split 1";
////        splitter->N3().Name_ = "Split 2";
//        system.add(splitter);
//
//        hypro::systemModule::ModulePtr split_isolator(new hypro::Friction(
//                "split_Isolator", &splitter->N2(), system.nodes_[3+nextNode].get()));
//        split_isolator->N2().Name_ = "Split Isolator end";
//        hypro::Friction *split_isolatorC = (hypro::Friction *) split_isolator.get();
//        split_isolatorC->L_ = isolatorLength;
//        split_isolatorC->D_ = 2*0.089;
//        split_isolatorC->Cf_ = isolatorFriction;
//        system.add(split_isolator);
//
//        hypro::systemModule::ModulePtr split_combustor(new hypro::EffComb<hypro::Friction>("split_combustor",system.nodes_[2 + nextNode].get(),system.nodes_[4 + nextNode].get()));
//        split_combustor->N2().Name_ = "Split Chamber End";
//        hypro::EffComb<hypro::Friction>* split_combustorC = (hypro::EffComb<hypro::Friction>*)split_combustor.get();
//        // Full GRI30
////        split_combustorC->Rcoeff_[0] = 2;
////        split_combustorC->Rcoeff_[3] = 1;
////        split_combustorC->Pcoeff_[5] = 2;
//        //        Reduced Species
//        split_combustorC->Rcoeff_[0] = 2;
//        split_combustorC->Rcoeff_[1] = 1;
//        split_combustorC->Pcoeff_[3] = 2;
//        split_combustorC->eta_ = combeff;
////        split_combustorC->chokingFeedback_ = injector;
//        split_combustorC->Cf_ = combustorFriction;
//        split_combustorC->D_ = 2*0.089;
//        split_combustorC->L_ = combustorLength;
//        system.add(split_combustor);
//
//        hypro::systemModule::ModulePtr mixer(new hypro::ParMixer("mixer",&split_isolator->N2(),system.nodes_[6].get(),&split_combustor->N2()));
//        system.add(mixer);

    } else if (combustionMode == 3){
        // ====================== Split EquilComb ======================
//        hypro::systemModule::ModulePtr splitter(new hypro::simpleSplitter("splitter",system.nodes_[0 + nextNode].get(),system.nodes_[1 + nextNode].get(),system.nodes_[2 + nextNode].get()));
//        splitter->N2().Name_ = "Split 1";
////        splitter->N3().Name_ = "Split 2";
//        system.add(splitter);
//
//        hypro::systemModule::ModulePtr split_isolator(new hypro::Friction(
//                "split_Isolator", &splitter->N2(), system.nodes_[3+nextNode].get()));
//        split_isolator->N2().Name_ = "Split Isolator end";
//        hypro::Friction *split_isolatorC = (hypro::Friction *) split_isolator.get();
//        split_isolatorC->L_ = isolatorLength;
//        split_isolatorC->D_ = 2*0.089;
//        split_isolatorC->Cf_ = isolatorFriction;
//        system.add(split_isolator);
//
//        hypro::systemModule::ModulePtr split_combustor(new hypro::EquilCombustion<hypro::Friction>("split_combustor",system.nodes_[2 + nextNode].get(),system.nodes_[4 + nextNode].get()));
//        split_combustor->N2().Name_ = "Split Chamber End";
//        hypro::EquilCombustion<hypro::Friction>* split_combustorC = (hypro::EquilCombustion<hypro::Friction>*)split_combustor.get();
//        split_combustorC->verbosity_ = 0;
//        split_combustorC->chokingFeedback_ = injector;
//        split_combustorC->Cf_ = combustorFriction;
//        split_combustorC->D_ = 2*0.089;
//        split_combustorC->L_ = combustorLength;
//        system.add(split_combustor);
//
////        hypro::systemModule::ModulePtr mixer(new hypro::ParMixer("mixer",&split_isolator->N2(),&system.N2(),&split_combustor->N2()));
//        hypro::systemModule::ModulePtr mixer(new hypro::ParMixer("mixer",&split_isolator->N2(),system.nodes_[6].get(),&split_combustor->N2()));
//        system.add(mixer);
    }



    system.calculate();
//    std::cout << system.nodes_[4+nextNode].get()->report() << std::endl;
//    std::cout << system.nodes_[3+nextNode].get()->report() << std::endl;

    std::cout << system.N1().report() << std::endl;
    std::cout << system.N2().report() << std::endl;
    hypro::Collection::OutType out = system.out();
    system.printOut(std::cout);
    vector<double> propellant = system.propellantMfr();

//    for (auto &massFlow : propellant){
//        std::cout << massFlow << ", ";
//    }

    std::cout << inlet.gasmodel_.get()->id() << ", " << system.N2().rho() << ", " << system.N2().getU() << ", " << system.N2().getTemp() << ", " << system.N2().getPress() << ", " << system.N2().mfr() << ", ";
    Cantera::compositionMap outlet = system.N2().getYMap();
    std::cout << outlet["H2"] << ", " << outlet["O2"] << ", " << outlet["N2"] << ", " << outlet["H2O"] << ", " << outlet["OH"] << ", " << outlet["O"] << ", " << outlet["H"] << ", " << outlet["NO"] << ", " << outlet["NO2"] << std::endl;

    return new hypro::systemModule(system);

}