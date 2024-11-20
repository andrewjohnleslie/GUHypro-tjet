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
//testSTATALTEXRamjet.cpp -- tests the STATALTEXRamjet specific impulses
#include "STATALTEXRamjet.h"
#include <gtest/gtest.h>
#include "US76.h"
#include "fmt/format.h"

///This test is to ensure that specific impulse predicted by the STATALTEXRamjet STX05 agree with
/**Figure 4.34 of <a href="http://oleg.lib.strath.ac.uk/R/?func=dbin-jump-full&object_id=26892">
 * "Toward automated design of Combined Cycle Propulsion." by Alessandro Mogavero</a>
 */
TEST(STATALTEXRamjet, STX05SpecificImpulse) {
double tolerance = 0.01;
    fmt::print("Here0");
    std::vector<double> mach;
    std::vector<double> altitude;
    std::vector<double> expectedSI = {1191.73, 1196.51, 1201.11, 1204.93, 1206.57, 1208.42, 1208.61, 1206.57, 1199.54,
                                      1191.73, 1183.17, 1171.49, 1148.56, 1123.96, 1096.07, 1066.46, 1032.15, 988.496,
                                      948.241, 929.256, 920.733, 910.943, 902.568, 894.764, 885.864, 881.626, 878.328,
                                      873.564, 870.564, 867.69, 864.986, 862.257, 859.7, 858.793, 857.877, 854.855,
                                      852.266, 850.359, 848.138, 847.245, 846.181, 845.183, 844.025, 843.575, 843.252,
                                      843.259, 843.7, 843.719, 843.452, 844.053, 845.588, 846.356, 847.115, 849.457,
                                      851.025, 852.725, 854.181, 854.694, 856.268, 857.413, 859.158, 861.073, 864.612,
                                      866.395, 868.896, 871.245, 873.64, 875.459, 877.884, 880.271, 882.847, 885.183,
                                      886.494, 887.797, 889.603, 890.788, 891.822, 893.561, 895.839, 897.065, 898.594,
                                      900.92, 902.695, 904.969, 907.295, 908.995, 910.508, 912.621, 914.238, 916.193,
                                      918.33, 919.519, 922.685, 924.533, 925.825, 926.875};

    std::string directory(QUOTE(STXRescource));
    std::ifstream ifs(directory + "/STX05TrajZM.csv");
    std::vector<std::vector<std::string>> data;

    fmt::print("Here1");
    while (ifs) {
        std::string s;
        if (!getline(ifs, s))
            break;

        istringstream ss(s);
        std::vector<std::string> record;

        while (ss) {
            std::string s;
            if (!getline(ss, s, ','))
                break;
            record.push_back(s);
        }
        data.push_back(record);
    }

    fmt::print("Here2");
    if (!ifs.eof())
        throw std::logic_error("Error: File STX05TrajZM.csv not found in search path: " + directory + "/");

    fmt::print("Here3");
    for (unsigned i = 0; i < data.size(); i++) {
        mach.push_back(stod(data[i][0]));
        altitude.push_back(stod(data[i][1]));
    }
    fmt::print("Here4");
    for (unsigned i = 0; i < data.size(); i++) {
        US76 atmo(altitude[i]);
        fmt::print("{},", i);
        double Thrust;
        std::vector<double> Mfr;
        hypro::Collection::OutType out;

        hypro::systemModule* sys = Ramjet(atmo.p(), atmo.T(), mach[i], 0, 5, Thrust, Mfr, out);
        Thrust = 0.96*sys->N2().A()*(sys->N2().rho()*std::pow(sys->N2().getU(),2) + sys->N2().getPress())
                 - sys->N1().A()*(sys->N1().rho()*std::pow(sys->N1().getU(),2) + sys->N1().getPress())
                 - sys->N1().getPress()*(sys->N2().A() - sys->N1().A());
        EXPECT_LT(1-std::fabs((Thrust / (9.81 * Mfr[0]))/expectedSI[i]), tolerance)<<"The expected specific thrust = " << Thrust / (9.81 * Mfr[0]);
        hypro::systemModule::deleteAll(sys);
    }
}

///This test is to ensure that thrust coefficient predicted by the STATALTEXRamjet STX05 agree with
/**Figure 4.33 of <a href="http://oleg.lib.strath.ac.uk/R/?func=dbin-jump-full&object_id=26892">
 * "Toward automated design of Combined Cycle Propulsion." by Alessandro Mogavero</a>
 */
TEST(STATALTEXRamjet, STX05ThrustCoefficient) {
    double tolerance = 0.01;

    std::vector<double> mach;
    std::vector<double> altitude;
    std::vector<double> expectedCT = {0.736754, 0.759074, 0.781632, 0.806774, 0.830456, 0.856015, 0.880768, 0.904689,
                                      0.922455, 0.940706, 0.959612, 0.978107, 0.989135, 0.993721, 0.994084, 0.994749,
                                      0.994898, 0.998081, 0.996076, 0.98544, 0.968521, 0.949431, 0.933371, 0.91861,
                                      0.90203, 0.894229, 0.888197, 0.879556, 0.874158, 0.868993, 0.864167, 0.859314,
                                      0.854787, 0.853216, 0.851601, 0.846372, 0.841965, 0.838709, 0.834989, 0.833507,
                                      0.83177, 0.830118, 0.828232, 0.827526, 0.827055, 0.82713, 0.827996, 0.828122,
                                      0.82774, 0.828853, 0.831614, 0.832997, 0.834403, 0.838579, 0.841395, 0.844452,
                                      0.847088, 0.848005, 0.850852, 0.852946, 0.856107, 0.859615, 0.866023, 0.869306,
                                      0.873905, 0.878227, 0.882675, 0.886067, 0.890588, 0.895068, 0.899918, 0.904461,
                                      0.907233, 0.909811, 0.913465, 0.915921, 0.918133, 0.921592, 0.926188, 0.928776,
                                      0.931975, 0.936716, 0.940327, 0.945017, 0.949825, 0.953365, 0.956567, 0.961027,
                                      0.964428, 0.968554, 0.97306, 0.975542, 0.982198, 0.986096, 0.988902, 0.991168};

    std::string directory(QUOTE(STXRescource));
    std::ifstream ifs(directory + "/STX05TrajZM.csv");
    std::vector<std::vector<std::string>> data;

    while (ifs) {
        std::string s;
        if (!getline(ifs, s))
            break;

        istringstream ss(s);
        std::vector<std::string> record;

        while (ss) {
            std::string s;
            if (!getline(ss, s, ','))
                break;
            record.push_back(s);
        }
        data.push_back(record);
    }
    if (!ifs.eof())
        throw std::logic_error("Error: File STX05TrajZM.csv not found in search path: " + directory + "/");

    for (unsigned i = 0; i < data.size(); i++) {
        mach.push_back(stod(data[i][0]));
        altitude.push_back(stod(data[i][1]));
    }

    for (unsigned i = 0; i < data.size(); i++) {
        US76 atmo(altitude[i]);

        double Thrust;
        std::vector<double> Mfr;
        hypro::Collection::OutType out;

        hypro::systemModule* sys = Ramjet(atmo.p(), atmo.T(), mach[i], 0, 5, Thrust, Mfr, out);
        Thrust = 0.96*sys->N2().A()*(sys->N2().rho()*std::pow(sys->N2().getU(),2) + sys->N2().getPress())
                 - sys->N1().A()*(sys->N1().rho()*std::pow(sys->N1().getU(),2) + sys->N1().getPress())
                 - sys->N1().getPress()*(sys->N2().A() - sys->N1().A());
        EXPECT_LT(1-std::fabs((Thrust/(0.5*sys->N1().rho()*std::pow(sys->N1().getU(), 2)*0.05432520556))/expectedCT[i]), tolerance)<<"The expected specific thrust = " << Thrust / (9.81 * Mfr[0]);
        hypro::systemModule::deleteAll(sys);
    }
}

///This test is to ensure that specific impulse predicted by the STATALTEXRamjet STX06 agree with
/**Figure 4.34 of <a href="http://oleg.lib.strath.ac.uk/R/?func=dbin-jump-full&object_id=26892">
 * "Toward automated design of Combined Cycle Propulsion." by Alessandro Mogavero</a>
 */
TEST(STATALTEXRamjet, STX06SpecificImpulse) {
    double tolerance = 0.01;

    std::vector<double> mach;
    std::vector<double> altitude;
    std::vector<double> expectedSI = {1186.69, 1183.24, 1170.43, 1169.02, 1165.94, 1165.12, 1165.11, 1168.92, 1170.63,
                                      1166.74, 1164, 1160.81, 1159.67, 1153.26, 1153.22, 1147.98, 1142.78, 1137.46,
                                      1127.24, 1098.97, 1068.94, 1030.24, 972.569, 930.263, 919.452, 909.973, 904.896,
                                      898.991, 893.193, 888.333, 883.287, 879.652, 875.338, 871.048, 865.06, 863.86,
                                      860.269, 857.018, 853.418, 848.987, 845.673, 841.165, 838.876, 836.448, 832.075,
                                      828.188, 824.236, 821.927, 816.067, 812.796, 809.527, 806.823, 804.136, 802.334,
                                      799.84, 797.58, 796.738, 794.469, 792.089, 790.656, 789.461, 788.165, 787.372,
                                      786.969, 786.047, 786.139, 786.22, 786.316, 786.703, 787.045, 787.263, 787.537,
                                      787.962, 788.809, 790.125, 792.119, 793.667};

    std::string directory(QUOTE(STXRescource));
    std::ifstream ifs(directory + "/STX06TrajZM.csv");
    std::vector<std::vector<std::string>> data;

    while (ifs) {
        std::string s;
        if (!getline(ifs, s))
            break;

        istringstream ss(s);
        std::vector<std::string> record;

        while (ss) {
            std::string s;
            if (!getline(ss, s, ','))
                break;
            record.push_back(s);
        }
        data.push_back(record);
    }
    if (!ifs.eof())
        throw std::logic_error("Error: File STX06TrajZM.csv not found in search path: " + directory + "/");

    for (unsigned i = 0; i < data.size(); i++) {
        mach.push_back(stod(data[i][0]));
        altitude.push_back(stod(data[i][1]));
    }

    for (unsigned i = 0; i < data.size(); i++) {
        US76 atmo(altitude[i]);

        double Thrust;
        std::vector<double> Mfr;
        hypro::Collection::OutType out;

        hypro::systemModule* sys = Ramjet(atmo.p(), atmo.T(), mach[i], 0, 6, Thrust, Mfr, out);
        Thrust = 0.96*sys->N2().A()*(sys->N2().rho()*std::pow(sys->N2().getU(),2) + sys->N2().getPress())
                 - sys->N1().A()*(sys->N1().rho()*std::pow(sys->N1().getU(),2) + sys->N1().getPress())
                 - sys->N1().getPress()*(sys->N2().A() - sys->N1().A());
        EXPECT_LT(1-std::fabs((Thrust / (9.81 * Mfr[0]))/expectedSI[i]), tolerance)<<"The expected specific thrust = " << Thrust / (9.81 * Mfr[0]);
        hypro::systemModule::deleteAll(sys);
    }
}

///This test is to ensure that thrust coefficient predicted by the STATALTEXRamjet STX06 agree with
/**Figure 4.33 of <a href="http://oleg.lib.strath.ac.uk/R/?func=dbin-jump-full&object_id=26892">
 * "Toward automated design of Combined Cycle Propulsion." by Alessandro Mogavero</a>
 */
TEST(STATALTEXRamjet, STX06ThrustCoefficient) {
    double tolerance = 0.01;

    std::vector<double> mach;
    std::vector<double> altitude;
    std::vector<double> expectedCT = {0.714594, 0.713861, 0.703825, 0.709052, 0.714447, 0.723952, 0.730497, 0.74785,
                                      0.762918, 0.771568, 0.78389, 0.7972, 0.811904, 0.826159, 0.842597, 0.859401,
                                      0.87795, 0.892393, 0.915511, 0.933819, 0.95028, 0.960332, 0.972071, 0.950759,
                                      0.929907, 0.912007, 0.902515, 0.891593, 0.881038, 0.872242, 0.863183, 0.85672,
                                      0.849125, 0.841605, 0.831186, 0.829129, 0.822946, 0.817423, 0.811304, 0.803837,
                                      0.798299, 0.790802, 0.787004, 0.78304, 0.775986, 0.769755, 0.763484, 0.759834,
                                      0.750649, 0.745557, 0.740496, 0.736347, 0.732242, 0.729488, 0.725683, 0.722289,
                                      0.721049, 0.717646, 0.71408, 0.711989, 0.710227, 0.708336, 0.707205, 0.706632,
                                      0.705314, 0.705516, 0.705733, 0.705934, 0.706596, 0.707189, 0.707563, 0.708042,
                                      0.708768, 0.7101, 0.712165, 0.715293, 0.717737};

    std::string directory(QUOTE(STXRescource));
    std::ifstream ifs(directory + "/STX06TrajZM.csv");
    std::vector<std::vector<std::string>> data;

    while (ifs) {
        std::string s;
        if (!getline(ifs, s))
            break;

        istringstream ss(s);
        std::vector<std::string> record;

        while (ss) {
            std::string s;
            if (!getline(ss, s, ','))
                break;
            record.push_back(s);
        }
        data.push_back(record);
    }
    if (!ifs.eof())
        throw std::logic_error("Error: File STX06TrajZM.csv not found in search path: " + directory + "/");

    for (unsigned i = 0; i < data.size(); i++) {
        mach.push_back(stod(data[i][0]));
        altitude.push_back(stod(data[i][1]));
    }

    for (unsigned i = 0; i < data.size(); i++) {
        US76 atmo(altitude[i]);

        double Thrust;
        std::vector<double> Mfr;
        hypro::Collection::OutType out;

        hypro::systemModule* sys = Ramjet(atmo.p(), atmo.T(), mach[i], 0, 6, Thrust, Mfr, out);
        Thrust = 0.96*sys->N2().A()*(sys->N2().rho()*std::pow(sys->N2().getU(),2) + sys->N2().getPress())
                 - sys->N1().A()*(sys->N1().rho()*std::pow(sys->N1().getU(),2) + sys->N1().getPress())
                 - sys->N1().getPress()*(sys->N2().A() - sys->N1().A());
        EXPECT_LT(1-std::fabs((Thrust/(0.5*sys->N1().rho()*std::pow(sys->N1().getU(), 2)*0.05432520556))/expectedCT[i]), tolerance)<<"The expected specific thrust = " << Thrust / (9.81 * Mfr[0]);
        hypro::systemModule::deleteAll(sys);
    }
}

///This test is to ensure that specific impulse predicted by the STATALTEXRamjet STX09 agree with
/**Figure 4.34 of <a href="http://oleg.lib.strath.ac.uk/R/?func=dbin-jump-full&object_id=26892">
 * "Toward automated design of Combined Cycle Propulsion." by Alessandro Mogavero</a>
 */
TEST(STATALTEXRamjet, STX09SpecificImpulse) {
    double tolerance = 0.01;

    std::vector<double> mach;
    std::vector<double> altitude;
    std::vector<double> expectedSI = {1043.22, 1033.21, 1030.45, 1017.71, 1012.39, 1009.05, 997.605, 989.276, 982.913,
                                      973.947, 966.435, 961.592, 956.843, 947.729, 940.031, 934.27, 932.048, 918.343,
                                      918.084, 919.306, 920.277, 921.213, 906.146, 908.758, 907.195, 903.52, 904.663,
                                      905.411, 904.248, 900.025, 897.095, 894.194, 887.483, 883.691};

    std::string directory(QUOTE(STXRescource));
    std::ifstream ifs(directory + "/STX09TrajZM.csv");
    std::vector<std::vector<std::string>> data;

    while (ifs) {
        std::string s;
        if (!getline(ifs, s))
            break;

        istringstream ss(s);
        std::vector<std::string> record;

        while (ss) {
            std::string s;
            if (!getline(ss, s, ','))
                break;
            record.push_back(s);
        }
        data.push_back(record);
    }
    if (!ifs.eof())
        throw std::logic_error("Error: File STX09TrajZM.csv not found in search path: " + directory + "/");

    for (unsigned i = 0; i < data.size(); i++) {
        mach.push_back(stod(data[i][0]));
        altitude.push_back(stod(data[i][1]));
    }

    for (unsigned i = 0; i < data.size(); i++) {
        US76 atmo(altitude[i]);

        double Thrust;
        std::vector<double> Mfr;
        hypro::Collection::OutType out;

        hypro::systemModule* sys = Ramjet(atmo.p(), atmo.T(), mach[i], 0, 9, Thrust, Mfr, out);
        Thrust = 0.96*sys->N2().A()*(sys->N2().rho()*std::pow(sys->N2().getU(),2) + sys->N2().getPress())
                 - sys->N1().A()*(sys->N1().rho()*std::pow(sys->N1().getU(),2) + sys->N1().getPress())
                 - sys->N1().getPress()*(sys->N2().A() - sys->N1().A());
        EXPECT_LT(1-std::fabs((Thrust / (9.81 * Mfr[0]))/expectedSI[i]), tolerance)<<"The expected specific thrust = " << Thrust / (9.81 * Mfr[0]);
        hypro::systemModule::deleteAll(sys);
    }
}

///This test is to ensure that thrust coefficient predicted by the STATALTEXRamjet STX09 agree with
/**Figure 4.33 of <a href="http://oleg.lib.strath.ac.uk/R/?func=dbin-jump-full&object_id=26892">
 * "Toward automated design of Combined Cycle Propulsion." by Alessandro Mogavero</a>
 */
TEST(STATALTEXRamjet, STX09ThrustCoefficient) {
    double tolerance = 0.01;

    std::vector<double> mach;
    std::vector<double> altitude;
    std::vector<double> expectedCT = {0.586179, 0.580005, 0.581278, 0.571522, 0.568872, 0.57007, 0.562187, 0.557515,
                                      0.553661, 0.547985, 0.544377, 0.54213, 0.539966, 0.534724, 0.530998, 0.530188,
                                      0.532954, 0.52368, 0.528363, 0.532339, 0.536763, 0.541331, 0.529506, 0.53814,
                                      0.543093, 0.544385, 0.55031, 0.556424, 0.559298, 0.563055, 0.56539, 0.567807,
                                      0.57247, 0.576104};

    std::string directory(QUOTE(STXRescource));
    std::ifstream ifs(directory + "/STX09TrajZM.csv");
    std::vector<std::vector<std::string>> data;

    while (ifs) {
        std::string s;
        if (!getline(ifs, s))
            break;

        istringstream ss(s);
        std::vector<std::string> record;

        while (ss) {
            std::string s;
            if (!getline(ss, s, ','))
                break;
            record.push_back(s);
        }
        data.push_back(record);
    }
    if (!ifs.eof())
        throw std::logic_error("Error: File STX09TrajZM.csv not found in search path: " + directory + "/");

    for (unsigned i = 0; i < data.size(); i++) {
        mach.push_back(stod(data[i][0]));
        altitude.push_back(stod(data[i][1]));
    }

    for (unsigned i = 0; i < data.size(); i++) {
        US76 atmo(altitude[i]);

        double Thrust;
        std::vector<double> Mfr;
        hypro::Collection::OutType out;

        hypro::systemModule* sys = Ramjet(atmo.p(), atmo.T(), mach[i], 0, 9, Thrust, Mfr, out);
        Thrust = 0.96*sys->N2().A()*(sys->N2().rho()*std::pow(sys->N2().getU(),2) + sys->N2().getPress())
                 - sys->N1().A()*(sys->N1().rho()*std::pow(sys->N1().getU(),2) + sys->N1().getPress())
                 - sys->N1().getPress()*(sys->N2().A() - sys->N1().A());
        EXPECT_LT(1-std::fabs((Thrust/(0.5*sys->N1().rho()*std::pow(sys->N1().getU(), 2)*0.05432520556))/expectedCT[i]), tolerance)<<"The expected specific thrust = " << Thrust / (9.81 * Mfr[0]);
        hypro::systemModule::deleteAll(sys);
    }
}
