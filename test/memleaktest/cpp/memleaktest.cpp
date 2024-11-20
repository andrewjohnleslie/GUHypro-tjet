
#include "Node.h"
//#include "thermoState.h"
//#include "cantera/IdealGasMix.h"
//#include "cantera/thermo/ThermoPhase.h"
//#include "cantera/thermo/ThermoFactory.h"
//#include "cantera/thermo/Phase.h"
//#include "thermoKineticState.h"
#include "systemModule.h"

#include "balanceMach.h"
#include "InjectionPlatePressure.h"
#include "EffIsenDuct.h"
#include "EffComb.h"
#include "Friction.h"
#include "AdaptedThroatInlet.h"
// #include "spdlog.h"
#include "fmt/format.h"


using namespace hypro;

//class thermoDumb
//{
//   public:
//       std::unique_ptr<Cantera::ThermoPhase> gasmodel_;
//       thermoDumb(const std::string &infile, std::string id_ = ""): gasmodel_(Cantera::newPhase(infile, id_)),
//                                                                   infile_(infile),
//                                                                   equilibrium(false)
//       {
//
//       }
//       ~thermoDumb(){}
//       std::string infile_;
//       bool equilibrium = false;
//
//};


 systemModule* Rocket() {


     bool rocket_verbosity = 1;

     // Create freeStream environment
     Node &freeStream = *new Node("/usr/share/cantera/data/gri30_highT.xml", "gri30");
     //Node &freeStream = *new Node("gri30_highT.xml", "gri30");
     freeStream.Name_ = "Freestream";
     if (rocket_verbosity != 1) {
         freeStream.warning_ = false;
     }
     Cantera::compositionMap initialComposition;
     initialComposition.emplace("O2", 0.05);
     initialComposition.emplace("H2", 1 - 0.05);
     freeStream.setTPX(300, 101324, initialComposition);
     freeStream.setU(300);

     // Create exit environment
     Node &exit = *new Node(freeStream);
     exit.Name_ = "Nozzle Exit";

     // Create system
     systemModule system("system", &freeStream, &exit);

     std::vector<double> empty_thermos(freeStream.X().size(), 0.0);

     vector<double> Xf = empty_thermos;
     const int numNodes = 4;
     static vector<double> A(numNodes);
     A = {1,0.9,0.9,0.9};

     freeStream.A(A[0]);

     system.initNodes(2, freeStream);
     for (int i = 0; i < 2; i++) {
         system.nodes_[i]->A(A[i + 1]);
     }
     exit.A(A[numNodes - 1]);


     systemModule::ModulePtr inlet(new EffIsenDuct("inlet",&system.N1(),system.nodes_[0].get()));
     inlet->N2().Name_ = "After Intake";
     EffIsenDuct* inletC = (EffIsenDuct*)inlet.get();
     inletC->etap_ = 1.0;
     inletC->etaT_ = 1.0;
     system.add(inlet);


     systemModule::ModulePtr nozzle2(new EffIsenDuct("convnozzle", &inlet->N2(), system.nodes_[1].get()));
     nozzle2->N2().Name_ = "Nozzle Throat";
     EffIsenDuct *nozzle2C = (EffIsenDuct *) nozzle2.get();
     nozzle2C->etap_ = 0.9;
     nozzle2C->etaT_ = 1;
     //nozzleC->choked_ = true;
     //propModule::chokingMach_ = 1.0;
     nozzle2C->verbosity_ = rocket_verbosity;
     system.add(nozzle2);


     hypro::systemModule::ModulePtr combustor(new hypro::EffComb<hypro::Friction>("combustor", system.nodes_[1].get(), &system.N1()));
     combustor->N2().Name_ = "Chamber End";
     hypro::EffComb<hypro::Friction> *combustorC =
             (hypro::EffComb<hypro::Friction> *)combustor.get();
     combustorC->Rcoeff_[0] = 2;
     combustorC->Rcoeff_[3] = 1;
     combustorC->Pcoeff_[5] = 2;
     combustorC->eta_ = 0.3;
     combustorC->verbosity_ = 0;
     // combustorC->chokingFeedback_ = injector;
     combustorC->Cf_ = 0.001264888478204313 + 0.0029981004934882564; //combustorFriction;
     combustorC->D_ = 2 * 0.093;										// 4*0.026*1/(2*1);
     combustorC->L_ = 0.126;
     system.add(combustor);

     system.verbosity_ = rocket_verbosity;
     system.calculate();

     return new systemModule(system);
 }



int main(int argc, char *argv[]) {

//     Node thisNode("/usr/share/cantera/data/gri30_highT.xml", "gri30");
//     Node newNode(thisNode);
//     thisNode.setTP(1000, 101325);
//     newNode.setTP(1000, 201325);

    // Cantera::IdealGasMix gas("/usr/share/cantera/data/gri30_highT.xml", "gri30");
    //  std::unique_ptr<Cantera::ThermoPhase> gasmodel_ (Cantera::newPhase("/usr/share/cantera/data/gri30_highT.xml", "gri30"));
//   thermoDumb thisstate("/usr/share/cantera/data/gri30_highT.xml", "gri30");
    // thermoState newstate(thisstate);
    // thermoDumb newdumb("/home/robert/repos/HyPro/test/memleaktest/gri30_highT.xml", "gri30");
//    thermoKineticState tks("/usr/share/cantera/data/gri30_highT.xml", "gri30");
//    thermoKineticState newtks(tks);
    // auto console = spdlog::stdout_color_mt("console");
    // console->set_level(spdlog::level::debug);
    // console->info("Welcome to spdlog");
     systemModule* sys = Rocket();
     systemModule::deleteAll(sys);


}