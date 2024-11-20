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
#include "Hyperion.h"

///Run Hyperion test case for all the test conditions
int main(int argc, char *argv[]) {

    int Machiterations = 9;
    int Altiterations = 4;
//    int Machiterations = 1;
//    int Altiterations = 1;

     double MachRange[9] = {2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0};
     double AltitudeRange[4] = {9144, 15240, 21336, 33528};

//    double MachRange[7] = {0.01, 0.5, 1.0, 1.5, 2.0, 2.5, 2.9};
//    double AltitudeRange[4] = {0, 6096, 12192, 18288};
//    double MachRange[1] = {0.01};
//    double AltitudeRange[1] = {0};

    for(int j=0; j<Altiterations; j++){
        double Altitude_ = AltitudeRange[j];
        US76 atmo(Altitude_);

        for(int i=0; i<Machiterations; i++){

            double Mach = MachRange[i];
            double Thrust;
            std::vector<double> Mfr;
            Collection::OutType out;
            fmt::print("Input: {}, {}, {}, {}\n", Altitude_, Mach, atmo.p(), atmo.T());

            systemModule* sys = Hyperion(atmo.p(),atmo.T(),Mach,1,1,Thrust,Mfr,out);
            systemModule::deleteAll(sys);

            std::cout  << "Thrust = " << Thrust << ", ISP = " << Thrust/(9.81*(Mfr[0])) <<  std::endl;

        }
    }
    return 0;
}
