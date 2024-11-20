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

#include "US76.h"
#include "STATALTEXRamjet.h"

int main(int argc, char *argv[]) {

    std::vector<double> mach;
    std::vector<double> altitude;

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
        Collection::OutType out;

        systemModule* sys = Ramjet(atmo.p(), atmo.T(), mach[i], 0, 5, Thrust, Mfr, out);

        Thrust = 0.96*sys->N2().A()*(sys->N2().rho()*std::pow(sys->N2().getU(),2) + sys->N2().getPress())
                 - sys->N1().A()*(sys->N1().rho()*std::pow(sys->N1().getU(),2) + sys->N1().getPress())
                 - sys->N1().getPress()*(sys->N2().A() - sys->N1().A());

        std::cout << "ALTITUDE: " << altitude[i] << ", Mach: " << mach[i] << std::endl;
        std::cout << "Thrust = " << Thrust << ", Specific Impulse = " << Thrust / (9.81 * Mfr[0]) << std::endl;
    }

    return 0;
}
