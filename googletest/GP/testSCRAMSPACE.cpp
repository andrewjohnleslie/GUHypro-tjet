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
//testCFD.cpp -- tests CFD case pressures and mach numbers
#include "SCRAMSPACE.h"
#include <gtest/gtest.h>
#include "fmt/format.h"

///The following test verifies the pressure output of the reacting Lorrain CFD case
/**(Figure 4.7 of <a href="http://oleg.lib.strath.ac.uk/R/?func=dbin-jump-full&object_id=26892">
 * "Toward automated design of Combined Cycle Propulsion." by Alessandro Mogavero</a>)
 */
TEST(SCRAMSPACE, PressureReacting) {
    double tolerance = 0.01;
    double p = 4100;
    double T = 370;
    double U = 2830;

    double Thrust;
    std::vector<double> Mfr;
    Collection::OutType out;

    hypro::systemModule* sys = SCRAMSPACE(p,T,U,0.4,0,Thrust,Mfr,out);
    EXPECT_LT(1-std::fabs(sys->N1().getPress()/4100), tolerance)<<"The pressure of inlet node = " << sys->N1().getPress();
    double internalPressure[4] = {32204.4, 42939.1, 51982.6, 111620};
    for (unsigned i = 0; i<sys->nodes_.size();i++) {
        EXPECT_LT(1-std::fabs(sys->nodes_[i].get()->getPress()/internalPressure[i]), tolerance)<<"The pressure of node " << i << " = " << sys->nodes_[i].get()->getPress();
    }
    EXPECT_LT(1-std::fabs(sys->N2().getPress()/19621.5), tolerance)<<"The pressure of exit node = " << sys->N2().getPress();
    hypro::systemModule::deleteAll(sys);
}

///The following test verifies the pressure output of the frozen Lorrain CFD case
/**(Figure 4.7 of <a href="http://oleg.lib.strath.ac.uk/R/?func=dbin-jump-full&object_id=26892">
 * "Toward automated design of Combined Cycle Propulsion." by Alessandro Mogavero</a>)
 */
TEST(SCRAMSPACE, PressureFrozen) {
    double tolerance = 0.01;
    double p = 4100;
    double T = 370;
    double U = 2830;

    double Thrust;
    std::vector<double> Mfr;
    Collection::OutType out;

    hypro::systemModule* sys = SCRAMSPACE(p,T,U,0,0,Thrust,Mfr,out);

    EXPECT_LT(1-std::fabs(sys->N1().getPress()/4100), tolerance)<<"The pressure of inlet node = " << sys->N1().getPress();
    double internalPressure[4] = {32204.4, 42939.1, 51982.6, 59170};
    for (unsigned i = 0; i<sys->nodes_.size();i++) {
        EXPECT_LT(1-std::fabs(sys->nodes_[i].get()->getPress()/internalPressure[i]), tolerance)<<"The pressure of node " << i << " = " << sys->nodes_[i].get()->getPress();
    }
    EXPECT_LT(1-std::fabs(sys->N2().getPress()/10786), tolerance)<<"The pressure of exit node = " << sys->N2().getPress();
}
