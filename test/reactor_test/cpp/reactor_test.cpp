#include "reactor_test.h"

int main(int argc, char *argv[]) {

    double p = 60000;
    //double p = 1;
    double T = 242.1912;
    double Mach = 4;

    double Thrust;
    vector<double> Mfr;
    hypro::Collection::OutType out;
    vector<double> emissions;


    Ramjet(p, T, Mach, Thrust, Mfr, out, emissions);

    double total_mass_flow = 0;
    vector<double> propellant_mfr;
    std::for_each(Mfr.begin(), Mfr.end(), [&] (double n) {
        total_mass_flow += n;
        if (n != 0){
            propellant_mfr.push_back(n);
        }
    });

    vector<double> emissions_mass = emissions*total_mass_flow;

    std::cout << std::endl << "Thrust = " << Thrust << ", Mfr = " << total_mass_flow << std::endl;
    std::cout << "Isp = " << Thrust/(9.81*total_mass_flow) << std::endl;
}