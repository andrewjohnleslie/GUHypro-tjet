#include "thermo/thermoKineticState.h"


void find_total_t(double T, double p, Cantera::compositionMap Y, double U){
    hypro::thermoKineticState &tks = *new hypro::thermoKineticState("/home/robert/repos/Cantera_HyPro/HyPro/test/scramjet_combustors/Burrows/gri30_highT.xml", "frc");
//    hypro::thermoKineticState &tks = *new hypro::thermoKineticState("/home/robert/repos/Cantera_HyPro/HyPro/test/scramjet_combustors/TotalT/nasa.xml", "gri30");
    tks.setTPY(T, p, Y);
    tks.setU(U);
//    std::cout  << tks.T0() << std::endl;
    std::cout  << ", " << tks.T0();
    std::cout << ", " << tks.p0();
//    std::cout << tks.report() << std::endl;
//
//    tks.equilibrium = true;
//    tks.setTP(T, p);
//
//    std::cout << tks.report() << std::endl;
//    std::cout << tks.T0() << std::endl;
}

