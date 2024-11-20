//
// Created by robert on 08/04/17.
//

#include "Node.h"

#include "core/systemModule.h"
#include "injectors/InjectionPlate.h"
#include "injectors/InjectionPlatePressure.h"
//#include "Combustion.h"
#include "solvers/EffIsenDuct.h"
#include "solvers/balanceMach.h"
#include "solvers/Friction.h"
#include "core/propModule.h"
#include "combustion/EquilCombustion.h"
#include "typedefs.h"
#include "returnData.h"


enum caseTag {
    SSME = 0,
    raptor = 1,
    raptor_vacuum = 2,
    merlin1D = 3,
    merlin1DVac=4,
    vulcain2 = 5,
    spaceliner_booster = 6,
    spaceliner_orbiter = 7,
    rs68 = 8,
    rs68a = 9,
    rl10b = 10,
    rl10a = 11,
};

using namespace hypro;
std::pair<bool, outputData>
Rocket(const double &p, const double &T, const double &Mach, const double throttle, caseTag engine, double &verbosity) {

    outputData returned;
    bool rocket_verbosity = 1;

    // Create freeStream environment
//    Node &freeStream = *new Node("/usr/share/cantera/data/gri30_highT.xml", "gri30");
//    Node &freeStream = *new Node("/home/robert/repos/HyPro/test/Rocket_test/gri30_highT_rocket.cti", "gri30");
//    Node &freeStream = *new Node("/home/robert/repos/HyPro/test/Rocket_test/JetSurf.cti", "gas");
    Node &freeStream = *new Node("/home/robert/repos/HyPro/test/Rocket_test/rocket_hydrocarbon.cti", "gri30");

    //Node &freeStream = *new Node("gri30_highT.xml", "gri30");
    freeStream.Name_ = "Freestream";
    if (rocket_verbosity != 1) {
        freeStream.warning_ = false;
    }
    Cantera::compositionMap initialComposition;
//    initialComposition.emplace("O2", 0.209);
//    initialComposition.emplace("N2", 1 - 0.209);
    initialComposition.emplace("C12H24", 1);
    freeStream.setTPX(T, p, initialComposition);
    vector<double> tempx = freeStream.X();
    // Hydrocarbon
    /*Node &freeStream = *new Node("/home/robert/repos/HyPro/test/Rocket_test/rp1.cti");
    //Node &freeStream = *new Node("/usr/share/cantera/data/gri30_highT.yaml", "gri30");
    freeStream.Name_ = "Freestream";
    if (rocket_verbosity != 1) {
        freeStream.warning_ = false;
    }
    std::vector<std::string> x = freeStream.XbyNames();

    Cantera::compositionMap initialComposition;
    initialComposition.emplace("O2", 0.209);
    initialComposition.emplace("N2", 1 - 0.209);
    freeStream.setTPX(T, p, initialComposition);
    std::cout << freeStream.report() << std::endl;*/

    std::vector<double> empty_thermos(freeStream.X().size(), 0.0);

    vector<double> Xf = empty_thermos;
    const int numNodes = 5;
    static vector<double> A(numNodes);
    //const int numNodes = 5;
    // double A[numNodes];

    double mfrMax;
    double injector_p;
    double injector_T;
    double o_f_ratio;

    switch (engine) {
        case SSME:
            // Space Shuttle Main Engine
            A = {0, 0.16, 0.16, 0.05, 4.15};
            mfrMax = 489;
            injector_p = 2e7;
            injector_T = 20;
            o_f_ratio = 6;
            Xf[0] = 1;
            Xf[3] = Xf[0] * o_f_ratio;
            std::cout << "Calculating SSME performance..." << std::endl;
            break;

        case raptor:
            // Raptor Engine
            // http://spaceflight101.com/spx/spacex-raptor/
            // Mixture Ratio: 3.8
            // Max Mass Flow Rate = 931.2 kg/s
            // Chamber Pressure = 30MPa
            // For SL:
            //      Area Ratio = 40
            //      Exit Diameter @ 1.7m
            // For Vac:
            //      Area Ratio = 200
            //      Exit Diameter = 4m
            // Predicted Thrust = 2.84988e+06
            // Reported Thrust = 3.050e+06
            A = {0, 0.16, 0.16, 0.056, 2.27};
            mfrMax = 931;
            injector_p = 3e7;
            injector_T = 20;
            o_f_ratio = 3.8;
            //Xf[13] = 1;
//            Xf[13] = 1;
//            Xf[3] = Xf[13] * o_f_ratio;
            Xf[6] = 1;
            Xf[144] = Xf[6] * o_f_ratio;
            std::cout << "Calculating Raptor Sea Level performance..." << std::endl;
            break;

        case raptor_vacuum:
            // Raptor Engine
            // http://spaceflight101.com/spx/spacex-raptor/
            // Mixture Ratio: 3.8
            // Max Mass Flow Rate = 931.2 kg/s
            // Chamber Pressure = 30MPa
            // For SL:
            //      Area Ratio = 40
            //      Exit Diameter @ 1.7m
            // For Vac:
            //      Area Ratio = 200
            //      Exit Diameter = 4m
            // Predicted Thrust = 3.24139e+06
            // Reported Thrust = 3.297e+06
            A = {0, 0.16, 0.16, 0.056, 11.2};
            mfrMax = 931;
            injector_p = 3e7;
            injector_T = 300;
            o_f_ratio = 3.8;
            Xf[13] = 1;
            Xf[3] = Xf[13] * o_f_ratio;
            std::cout << "Calculating Raptor Vacuum performance..." << std::endl;
            break;

        case merlin1D:
            // RP1 http://pubs.acs.org/doi/pdf/10.1021/ef900216z
            // alpha-methyldecalin 0.354 C11H20
            // 5-methylnonane 0.150 C10H20
            // n-dodecane 0.183 C12H26
            // heptylcyclohexane 0.313 C13H26
            //
//            A = {0, 0.16, 0.16, 0.042, 0.9};
            A = {0, 0.196349, 0.196349, 0.049087385, 0.78539816};
            mfrMax = 300;
            injector_p = 11.0e6;
            injector_T = 215;
            o_f_ratio = 2.36;
            //Xf[13] = 1;
            //Xf[3] = Xf[13] * o_f_ratio;
            Xf[53] = 1/(1 + o_f_ratio);
            Xf[3] = o_f_ratio / (1 + o_f_ratio);
            std::cout << "Calculating Merlin1D Vacuum performance..." << std::endl;
            break;

        case merlin1DVac:
            // RP1 http://pubs.acs.org/doi/pdf/10.1021/ef900216z
            // alpha-methyldecalin 0.354 C11H20
            // 5-methylnonane 0.150 C10H20
            // n-dodecane 0.183 C12H26
            // heptylcyclohexane 0.313 C13H26
            //
//            A = {0, 0.16, 0.16, 0.042, 0.9};
            A = {0, 0.196189605223361, 0.196189605223361, 0.04904740130584, 8.09282121546365};
            mfrMax = 292;
            injector_p = 11.0e6;
            injector_T = 215;
            o_f_ratio = 2.36;
            //Xf[13] = 1;
            //Xf[3] = Xf[13] * o_f_ratio;
            Xf[53] = 1/(1 + o_f_ratio);
            Xf[3] = o_f_ratio / (1 + o_f_ratio);
            std::cout << "Calculating Merlin1D Vacuum performance..." << std::endl;
            break;

        case vulcain2:
            // http://www.astronautix.com/v/vulcain2.html
            // Diameter = 2.10, exit area = 3.46m
            // Area Ratio = 61.5
            A = {0, 0.16, 0.16, 0.062, 3.46};
            mfrMax = 320;
            //mfrMax = 500;
   //         injector_p = 1.20e7;
            injector_p = 1.150e7;
            injector_T = 20;
            o_f_ratio = 6.1;
            Xf[0] = 1;
            Xf[3] = Xf[0] * o_f_ratio;
            std::cout << "Calculating Vulcain 2 performance..." << std::endl;
            break;

        case spaceliner_booster:
            // "Staged Combustion Cycle Rocket Engine Design Trade-Offs for Future Advanced Passenger Transport" by Sippel, Yamashiro, Cremaschi, 2012
            injector_p = 1.6e7;
            injector_T = 20;
            o_f_ratio = 6;
            mfrMax = 384.5;
            break;


        case spaceliner_orbiter:
            injector_p = 1.6e7;
            injector_T = 20;
            o_f_ratio = 6;
            mfrMax = 384.5;
            break;

        case rs68:
            // http://www.spaceflight101.net/delta-iv-heavy-rs-68a.html
            // https://pdfs.semanticscholar.org/bc61/88db410e8ca34f4cb6017b0f86f53d744ed1.pdf
            // Engine Diameter = 2.43m
            // Exit Area = 4.63
            // Nozzle Ratio = 21.5
            A = {0, 0.6, 0.6, 0.216, 4.63};
           // A = {0, 0.6, 0.6, 0.112, 4.63};
            injector_p = 1.01e+7;
            //injector_p = 2.00e+7;
            injector_T = 20;
            o_f_ratio = 6;
            mfrMax = 900;
            Xf[0] = 1;
            Xf[3] = Xf[0] * o_f_ratio;
            //std::cout << "Calculating RS-68 performance..." << std::endl;
            break;

        case rs68a:
            // http://www.spaceflight101.net/delta-iv-heavy-rs-68a.html
            // Engine Diameter = 2.43m
            // Exit Area = 4.63
            // Nozzle Ratio = 21.5
            A = {0, 0.8628274, 0.8628274, 0.215706866, 4.637697615};
//            // A = {0, 0.6, 0.6, 0.112, 4.63};
//            injector_p = 1.0559e+7;
            injector_p = 9.7e6;
//            injector_p = 19.6e6;
            injector_T = 20;
            o_f_ratio = 5.97;
//            mfrMax = 893.66; //FROZEN
//                    mfrMax = 882.22; //EQUIL
                     mfrMax = 890;
//            mfrMax = 650;
            Xf[0] = 1/(1 + o_f_ratio);
            Xf[3] = o_f_ratio / (1 + o_f_ratio);
            //std::cout << "Calculating RS-68 performance..." << std::endl;
            break;

        case rl10b:
            // http://spaceflight101.com/spacerockets/delta-iv-medium-42/
            //https://www.rocket.com/files/aerojet/documents/Capabilities/PDFs/RL10%20data%20sheet%20Feb%202016.pdf
            // De = 2.1463
            // Ae =  3.618
            // A* = 0.0144721
            // m_dot = 24.088
            // Area Ratio = 250
            A = {0,100,100,0.0134595,3.836}; //IAC
//            A = {0,0.05383808,0.05383808,0.01345952,3.83596317};
//            A = {0,0.05119674,0.05119674,0.012799185,3.647768}; //REAL?
//            injector_p = 4.6e6 ;
            injector_p = 4.413e6 ;
//            injector_p = 4.413e6 ;
            injector_T = 20 ;
            o_f_ratio = 5.88;
//            mfrMax = 30;
//        mfrMax = 24;
        mfrMax = 27;
//            Xf[0] = 1;
//            Xf[3] = Xf[0] * o_f_ratio;
            Xf[0] = 1/(1 + o_f_ratio);
            Xf[3] = o_f_ratio / (1 + o_f_ratio);
            break;

        case rl10a:
            A = {0,0.05119674,0.05119674,0.012799185, 1.075131546};
            injector_p = 4.206e6 ;
            injector_T = 30 ;
            o_f_ratio = 5.5;
            mfrMax = 23;
            Xf[0] = 1/(1 + o_f_ratio);
            Xf[3] = o_f_ratio / (1 + o_f_ratio);
    }

    // Create exit environment
    Node &exit = *new Node(freeStream);
    exit.Name_ = "Nozzle Exit";

    // Create system
    systemModule system("system", &freeStream, &exit);

    freeStream.A(A[0]);
    system.initNodes(3, freeStream);
    for (int i = 0; i < 3; i++) {
        system.nodes_[i]->A(A[i + 1]);
    }
    exit.A(A[numNodes - 1]);

    systemModule::ModulePtr injector(new InjectionPlatePressure("injector", &system.N1(), system.nodes_[0].get()));
    injector->N2().Name_ = "Injector Exit";
    injector->verbosity_ = rocket_verbosity;
    InjectionPlatePressure *injectorC = (InjectionPlatePressure *) injector.get();
    injectorC->mfrMax_ = mfrMax;
    injectorC->p_ = injector_p;
    //vector<double> Xf = empty_thermos;
    //double o_f_ratio = 6;
    //Xf[0] = 1;
    //Xf[3] = Xf[0]*o_f_ratio;
    injectorC->Y_ = Xf;
    injectorC->T_ = injector_T;
    injectorC->unreduce();
    injectorC->verbosity_ = rocket_verbosity;
    system.add(injector);

    systemModule::ModulePtr combustor(
            new EquilCombustion<balanceMach>("combustor", system.nodes_[0].get(), system.nodes_[1].get()));
    combustor->N2().Name_ = "Combustor Outlet";
    EquilCombustion<balanceMach> *combustorC = (EquilCombustion<balanceMach> *) combustor.get();
    combustorC->verbosity_ = rocket_verbosity;
    combustorC->chokingFeedback_ = injector;
    combustorC->solver_verbosity = false;
    system.add(combustor);

    systemModule::ModulePtr nozzle(new EffIsenDuct("convnozzle", system.nodes_[1].get(), system.nodes_[2].get()));
    nozzle->N2().Name_ = "Nozzle Throat";
    EffIsenDuct *nozzleC = (EffIsenDuct *) nozzle.get();
    nozzleC->etap_ = 1;
    nozzleC->etaT_ = 1;
    nozzleC->chokingFeedback_ = injector;
    nozzleC->N2().setEquilibrium(true);
    nozzleC->verbosity_ = false;
    system.add(nozzle);

    systemModule::ModulePtr divnozzle(new EffIsenDuct("divnozzle", system.nodes_[2].get(), &system.N2()));
    divnozzle->N2().Name_ = "Nozzle Exit";
    EffIsenDuct *divnozzleC = (EffIsenDuct *) divnozzle.get();
    divnozzleC->etap_ = 1;
//    propModule::chokingMach_ = 0.9;
//    propModule::chokingMach_ = 0.9;
    propModule::chokingMach_ = 0.9;
    divnozzleC->etaT_ = 1;
    divnozzleC->choked_ = true;
    divnozzleC->N2().setEquilibrium(true);
    divnozzleC->chokingFeedback_ = injector;
    divnozzleC->verbosity_ = false;
    system.add(divnozzle);

    system.verbosity_ = rocket_verbosity;
    system.calculate();
    system.printOut(std::cout);
//    std::cout << system.N2().report() << endl;

    returned.thrust = system.thrust();
    returned.Mfr = system.propellantMfr();
    returned.out = system.out();
    //emissions = system.N2().X();

    returned.total_flow = system.N2().rho()*system.N2().getU()*system.N2().A();
    //std::cout << system.N2().report() << std::endl;
//    emissions.clear();
    fmt::print(system.N2().report());
    returned.emissions.push_back(system.N2().Y("H2"));
    returned.emissions.push_back(system.N2().Y("H"));
    returned.emissions.push_back(system.N2().Y("O"));
    returned.emissions.push_back(system.N2().Y("O2"));
    returned.emissions.push_back(system.N2().Y("OH"));
    returned.emissions.push_back(system.N2().Y("H2O"));
    returned.emissions.push_back(system.N2().Y("HO2"));
    returned.emissions.push_back(system.N2().Y("H2O2"));

    returned.sys = new systemModule(system);
//    return new systemModule(system);
    return {true, returned};
}
