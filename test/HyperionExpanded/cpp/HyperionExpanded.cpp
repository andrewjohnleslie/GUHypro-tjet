/*!
 *  \author     Alessandro Mogavero
 *  \n
 *  \email      alessandro.mogavero@strath.ac.uk
 *  \version    1.0
 *  \copyright  Copyright 2016 Alessandro Mogavero
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
* along with HyPro.  If not, see <http://www.gnu.org/licenses/>. */

#include "HyperionExpanded.h"
#include "US76.h"
//#include "fenv.h"
#include <boost/math/tools/roots.hpp>


void find_dynamic_pressure() {
    /*
        Routine finds the altitude at which the input mach number and the dynamic pressure are true.
    */
    double Dynamic_pressure = 95760.5; //2000psf
    double mach = 10;
    double H2;

    double Hmin = 0;
    double Hmax = 80000;

//    double factor = 1.25;

    const boost::uintmax_t maxit = 30;
    boost::uintmax_t it = maxit;
//    bool is_rising = true;
    int digits = std::numeric_limits<double>::digits;
    int get_digits = digits - 3;
    boost::math::tools::eps_tolerance<double> tol(get_digits);

    auto G_err = [mach, Dynamic_pressure](double H) {
        US76 atmo(H);
        Node &freeStream = *new Node("/usr/share/cantera/data/gri30_highT.yaml", "gri30");
        Cantera::compositionMap initialComposition;
        initialComposition.emplace("O2", 0.209);
        initialComposition.emplace("N2", 1 - 0.209);

        freeStream.setTPX(atmo.T(), atmo.p(), initialComposition);


        double qcalc = (freeStream.rho() * mach * freeStream.geta() * mach * freeStream.geta()) / 2;
        std::cout << "Altitude: " << H << std::endl;
        std::cout << "Density: " << freeStream.rho() << std::endl;
        std::cout << "Speed of Sound: " << freeStream.geta() << std::endl;
        std::cout << "Dynamic Pressure: " << qcalc << std::endl;

        delete &freeStream;

        std::cout << "Diff: " << qcalc - Dynamic_pressure << std::endl;
        return qcalc - Dynamic_pressure;
    };

    //std::pair<double, double> res = boost::math::tools::bracket_and_solve_root(G_err, guess, factor, is_rising, tol,
    //                                                                           it);
    std::pair<double, double> res = boost::math::tools::toms748_solve(G_err, Hmin, Hmax, tol, it);

    H2 = res.first + (res.second - res.first) / 2;

    std::cout << H2 << std::endl;

    US76 atmo(H2);
    Node &freeStream = *new Node("/usr/share/cantera/data/gri30_highT.yaml", "gri30");
    Cantera::compositionMap initialComposition;
    initialComposition.emplace("O2", 0.209);
    initialComposition.emplace("N2", 1 - 0.209);

    freeStream.setTPX(atmo.T(), atmo.p(), initialComposition);


    double qcalc = (freeStream.rho() * mach * freeStream.geta() * mach * freeStream.geta()) / 2;
    std::cout << "Density: " << freeStream.rho() << std::endl;
    std::cout << "Speed of Sound: " << freeStream.geta() << std::endl;
    std::cout << "Dynamic Pressure: " << qcalc << std::endl;


    std::cout << atmo.T() << std::endl;
    std::cout << atmo.p() << std::endl;
    std::cout << atmo.rho() << std::endl;

    //std::cout << N2_->report() << std::endl;
    // std::cout << "M2: " << M2 << std::endl;
}

void run_model(const int verbosity, const int mode, const bool labels, const bool printEmissions,
               const bool printMassFlowRates, const vector<double> AltitudeRange, const vector<double> MachRange,
               const vector<double> throttleRange) {

//  ##### Step 4 - Run Hyperion through input data provided #####
    for (auto &Throttle : throttleRange) {
        for (auto &Altitude : AltitudeRange) {
            US76 atmo(Altitude);
            for (auto &Mach : MachRange) {

                fmt::print("Input: {}, {}, {}, {}, {}, {}\n", atmo.p(), atmo.T(), Mach, Throttle, mode, verbosity);
                auto[ok, data] = Hyperion(atmo.p(), atmo.T(), Mach, Throttle, mode, verbosity);
                //
//                if (Mach == 0 &&  Altitude == 0) {
//                //
//                //     // create and open geta character archive for output
//                //     std::ofstream ofs("Hyperion.xml");
//                //     sys->save(ofs);
//                //     ofs.close();
                data.sys->saveJson("Hyperion.hyp");
//
//			std::shared_ptr<propModule> sys1 = propModule::loadJson("Hyperion.hyp");
//
//			sys1->unreduce();
//			sys1->calculate();
//			((systemModule*)sys1.get())->printOut(std::cout);
//
//			std::stringstream ss1;
//			sys1->saveJson(ss1);
//			std::cout << ss1.str() << std::endl;
//                 }
                systemModule::deleteAll(data.sys);

                if (labels == true) {
                    fmt::print("Altitude = {}, Mach = {}, Pressure = {}, Temperature = {}, Thrust = {}, ISP = {}, "
                               "Ct = {}, Mfr = {}, q = {}\n", Altitude, Mach, atmo.p(), atmo.T(), data.thrust,
                               data.thrust / (9.81 * (data.Mfr[0] + data.Mfr[3])),
                               data.thrust / (data.dynamic_pressure * 2.51), data.Mfr[0] + data.Mfr[3],
                               data.dynamic_pressure);

                } else {
                    if (mode == 0) {
                        fmt::print("{}, {}, {}, {}, {}, {}, {}, {}\n", Altitude, Mach, atmo.p(), atmo.T(),
                                   data.thrust, data.thrust / (9.81 * (data.Mfr[0] + data.Mfr[3])), "Nan",
                                   data.Mfr[0] + data.Mfr[3]);
                    } else {
                        fmt::print("{}, {}, {}, {}, {}, {}, {}, {}\n", Altitude, Mach, atmo.p(), atmo.T(),
                                   data.thrust, data.thrust / (9.81 * (data.Mfr[0] + data.Mfr[3])),
                                   data.thrust / (data.dynamic_pressure * 2.51), data.Mfr[0] + data.Mfr[3]);
                    }
                }

                if (printEmissions == true) {
                    vector<double> emissions_mass = data.emissions * data.total_flow;
                    for (const double &emission : data.emissions * data.total_flow) { fmt::print("{}, ", emission); }
                    fmt::print("\n");
                }

                if (printMassFlowRates == true) {
                    for (const double &massFlow : data.Mfr) { fmt::print("{}, ", massFlow); }
                    fmt::print("\n");
                }

            }
        }
    }

    /*
    Node &freeStream = *new Node("/usr/share/cantera/data/gri30_highT.yaml", "gri30");
    freeStream.Name_ = "Freestream";

    Cantera::compositionMap initialComposition;
    initialComposition.emplace("O2", 0.209);
    initialComposition.emplace("N2", 1 - 0.209);

    for (int t = 0; t < 100; t++) {
        US76 atmo(t*1000);

        freeStream.setTPX(atmo.T(),atmo.p(),initialComposition);
        freeStream.M();
    }*/

    //find_dynamic_pressure();
}

///Run Hyperion test case for all the test conditions
int main(int argc, char *argv[]) {

//    ##### Step 1 - Select Engine Mode, Mach Range and Altitude Range #####
//
//    Nominal test values:
//    Ejector:
//    SCCREAM (1997)
//    mode = 0
//    std::vector<double> MachRange = {0, 0.5, 1, 1.5, 2, 2.5, 2.9};
//    std::vector<double> AltitudeRange = {0, 6096, 12192, 18288};

//    Ramjet:
//    SCCREAM (1997)
//    mode = 1;
    std::vector<double> MachRange = {2.0, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6};
    std::vector<double> AltitudeRange = {9144, 15240, 21336, 33528}; // AIAA-1997
//    std::vector<double> AltitudeRange = {9144, 15240, 21336, 27432}; // AIAA-2001
//    std::vector<double> AltitudeRange = {0,500,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000,6500,7000,7500,8000,8500,9000,9500,10000,10500,11500,12000,12500,13000,13500,14000,14500,15000,15500, 16000,16500,17000,17500,18000,18500,19000,19500,20000};

//    Scramjet:
//    mode = 2;
//    std::vector<double> MachRange =  {6.0, 7.0, 8.0, 9.0, 10.0};
//    std::vector<double> AltitudeRange =  {15240, 21336, 27432};

//    Rocket
//    mode = 3;
//    std::vector<double> MachRange = {0}
//    std::vector<double> AltitudeRange = {0,10000,20000,30000,40000,50000,60000,70000,80000,90000,100000)

//    SINGLE:
//    std::vector<double> MachRange = {3.5};
//    std::vector<double> AltitudeRange = {17373.6};


//  ##### Step 2 - Select range of throttles to calculate over #####
    std::vector<double> throttleRange = {1};

//  ##### Step 3 - Input mode select, verbosity setting and printing settings
    run_model(1, 0, true, true, true, AltitudeRange, MachRange, throttleRange);

    return 0;
}
