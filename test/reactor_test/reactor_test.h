//
// Created by robert on 08/04/17.
//

#include "core/Node.h"
#include "core/systemModule.h"
#include "injectors/InjectionPlate.h"
#include "injectors/InjectionPlatePressure.h"
//#include "Combustion.h"
#include "solvers/EffIsenDuct.h"
#include "solvers/balanceMach.h"
#include "solvers/Friction.h"
#include "core/propModule.h"
#include "combustion/EquilCombustion.h"
#include "combustion/CombustionReactor.h"
#include "typedefs.h"
#include "Wedge.h"
#include "inletModels/AdaptedThroatInlet.h"
#include "injectors/InjectionPhi.h"
#include "combustion/EffComb.h"
#include "AdaptedNozzle.h"

//using namespace hypro;

vector<double> CreateEmptyX(hypro::Node inputX){
    vector<double> empty_X = inputX.X();
    empty_X[3] = 0;
    empty_X[47] = 0;
    return empty_X;
}


hypro::systemModule *Ramjet(const double &p, const double &T, const double &Mach, double &Thrust, vector<double> &Mfr,
       hypro::Collection::OutType &out, vector<double> &emissions) {

    bool rocket_verbosity = 0;
    bool reactor_on = true;
    double foot = 0.3048; //conversion feet to meters
    double lbm = 0.45359; //conversion libre to kg


    // Create freeStream environment
    hypro::Node &freeStream = *new hypro::Node("/usr/share/cantera/data/gri30_highT.yaml", "gri30");
    freeStream.Name_ = "Freestream";
    if (rocket_verbosity != 1) {
        freeStream.warning_ = false;
    }
    Cantera::compositionMap initialComposition;
    initialComposition.emplace("O2", 0.209);
    initialComposition.emplace("N2", 1 - 0.209);
    freeStream.setTPX(T, p, initialComposition);
    freeStream.setU(Mach * freeStream.geta());
    freeStream.Name_ = "Free Stream";
    std::cout << freeStream.report() << std::endl;

    // Create empty X vector (useful for fuels etc.)
    std::vector<double> empty_thermos = CreateEmptyX(freeStream);

    // Create Modified freestream
    hypro::Node &freeStreamMod = *new hypro::Node(freeStream);
    freeStreamMod.Name_ = "Pre Intake";

    hypro::Wedge forebody("forebody", &freeStream, &freeStreamMod);
    forebody.delta_ = 8.0;
    forebody.calculate();

    // Create exit environment
    hypro::Node &exit = *new hypro::Node(freeStream);
    exit.Name_ = "Nozzle Exit";

    // Create system
    hypro::systemModule system("system", &freeStream, &exit);

    //Nodes                 N1_    I0    I1       0     1       2     3      4        5      6     7     N2_
    const int numNodes = 12;
    double A[numNodes] = {27.0, 27.0, 0.25 * 27, 8.24, 11.25, 11.25, 0.0, 11.25 - 8.24, 22.5, 22.5, 22.5, 44.3};
    /*name = {'Pre Intake','Intake','Throat','Pinch Point','Primary','Mixer End',...
        'Dummy','Rocket Outlet','End Diffuser','Injection','End Chamber','Nozzle'};*/
    double foot2 = std::pow(foot, 2.0);
    freeStreamMod.A(A[0] * foot2);
    system.initNodes(8, freeStreamMod);
    for (int i = 0; i < 8; i++) {
        system.nodes_[i]->A(A[i + 3] * foot2);
    }
    exit.A(A[numNodes - 1] * foot2);
    exit.Amax_ = A[numNodes - 1] * foot2;

    // ============================ CREATE MODULES ==============================
    // ################################ Inlet
    hypro::systemModule::ModulePtr inlet(new hypro::AdaptedThroatInlet("Inlet", &system.N1(), system.nodes_[0].get()));
    inlet->N2().Name_ = "Pinch Point";
    hypro::AdaptedThroatInlet *inletC = (hypro::AdaptedThroatInlet *) inlet.get();
    inletC->nodes_[0]->A(A[1] * foot2);
    inletC->nodes_[0]->Amax_ = inletC->nodes_[0]->A();
    inletC->nodes_[0]->Name_ = "Intake";
    inletC->nodes_[1]->A(1.055 * A[2] * foot2);
    inletC->nodes_[1]->Amax_ = inletC->nodes_[1]->A();
    inletC->nodes_[1]->Name_ = "Throat";
    inlet->verbosity_ = 0;
    inlet->chokingFeedback_ = inlet;
    system.add(inlet);

    // ################################ Primary / ejector
    hypro::systemModule::ModulePtr postPinch(new hypro::isentropicDuct("Post pinch point", system.nodes_[0].get(), system.nodes_[1].get()));
    system.add(postPinch);

    hypro::systemModule::ModulePtr primary(new hypro::EffIsenDuct("primary", system.nodes_[1].get(), system.nodes_[2].get()));
    primary->chokingFeedback_ = inlet;
    hypro::EffIsenDuct *primaryC = (hypro::EffIsenDuct *) primary.get();
    primaryC->etap_ = 0.725;
    system.add(primary);

    // ################################ Post Combustor
    hypro::systemModule::ModulePtr diffuser(new hypro::EffIsenDuct("diffuser", system.nodes_[2].get(), system.nodes_[5].get()));
    diffuser->N2().Name_ = "Diffuser End";
    hypro::EffIsenDuct *diffuserC = (hypro::EffIsenDuct *) diffuser.get();
    diffuserC->etap_ = 0.92;
    diffuserC->etaT_ = 1.0;
    system.add(diffuser);

    hypro::systemModule::ModulePtr secondary(
            new hypro::InjectionPhi<hypro::balanceMach>("secondary", system.nodes_[5].get(), system.nodes_[6].get()));
    secondary->N2().Name_ = "Injection";
    hypro::InjectionPhi<hypro::balanceMach> *secondaryC = (hypro::InjectionPhi<hypro::balanceMach> *) secondary.get();
    secondaryC->PhiMax_ = 1; //TODO Add throttle for injector here
    secondaryC->Phi_ = 1;
    //secondaryC->mfrMax_ = *throttle;
    vector<double> Xf = empty_thermos;
    Xf[0] = 1;
    secondaryC->Xf_ = Xf;
    secondaryC->Tf_ = 20;

    secondaryC->Rcoeff_[0] = 2;
    secondaryC->Rcoeff_[3] = 1;
    secondary->chokingFeedback_ = secondary;
    system.add(secondary);

    // ############################### Combustor
    // Reactors
    if (reactor_on == true) {
        hypro::systemModule::ModulePtr combustor(
                new hypro::CombustionReactor("combustor", system.nodes_[6].get(), system.nodes_[7].get()));
    combustor->N2().Name_ = "Chamber End";
    //combustor->chokingFeedback_ = secondary;
    hypro::CombustionReactor *combustorC = (hypro::CombustionReactor *) combustor.get();
//    combustorC->Rcoeff_[0] = 2;
//    combustorC->Rcoeff_[3] = 1;
//    combustorC->Pcoeff_[5] = 2;
//    combustorC->eta_ = 0.8;
//    combustorC->Cf_ = 0.01;
//    combustorC->D_ = 1.631;
    combustorC->length = 10.0;
    system.add(combustor);
    }else{

    // Equilibrium
        hypro::systemModule::ModulePtr combustor(
        new hypro::EquilCombustion<hypro::balanceMach>("combustor", system.nodes_[6].get(), system.nodes_[7].get()));
    combustor->N2().Name_ = "Combustor Outlet";
        hypro::EquilCombustion<hypro::balanceMach> *combustorC = (hypro::EquilCombustion<hypro::balanceMach> *) combustor.get();
    combustorC->verbosity_ = rocket_verbosity;
    combustorC->chokingFeedback_ = inlet;
    system.add(combustor);
    }


    // ################################ Nozzle
    hypro::systemModule::ModulePtr nozzle(new hypro::AdaptedNozzle<hypro::EffIsenDuct>("nozzle",system.nodes_[7].get(),&system.N2()));
    hypro::AdaptedNozzle<hypro::EffIsenDuct>* nozzleC = (hypro::AdaptedNozzle<hypro::EffIsenDuct>*)nozzle.get();
    nozzleC->choked_ = true;
    /* in case of inlet spilling (especially in supersonic) better to adapt the
        nozzle to the inlet conditions.*/
    nozzleC->freeStream_ = &inlet->N1();
    nozzleC->etap_ = 1.0;
    system.add(nozzle);


    // ============================ CALCULATE ==============================
    system.verbosity_ = rocket_verbosity;
    system.calculate();
    system.printOut(std::cout);
    std::cout << system.N2().report() << endl;

    Thrust = system.thrust();
    Mfr = system.propellantMfr();
    out = system.out();
    emissions = system.N2().X();

    return new hypro::systemModule(system);
}

