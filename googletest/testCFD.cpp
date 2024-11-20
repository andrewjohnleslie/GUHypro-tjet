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
//#include "CFD.h"
#include <gtest/gtest.h>

///This test is to ensure that the provided inlet conditions are the same as those tested in
/**<a href="http://oleg.lib.strath.ac.uk/R/?func=dbin-jump-full&object_id=26892">
 * "Toward automated design of Combined Cycle Propulsion." by Alessandro Mogavero</a>
 */
/*TEST(CFD, inletConds) {
    double tolerance = 0.01;
    std::unique_ptr<hypro::systemModule> sys = CFD(CASE1tuned,1.0);
    EXPECT_LT(1-std::fabs(sys->N1().gamma()/1.3735), tolerance)<<"The inlet gamma = " << sys->N1().gamma();
    EXPECT_LT(1-std::fabs((sys->N1().getT0())/620.125), tolerance)<<"The inlet total temp = " << sys->N1().T0();
    EXPECT_LT(1-std::fabs((sys->N1().mfr()*2/360)/4.13019), tolerance)<<"The inlet flow rate = " << (sys->N1().mfr()*2/360);
    EXPECT_LT(1-std::fabs((1.5*std::pow(sys->N1().getU()*0.1,2.0))/314.32), tolerance)<<"The inlet k = " << (1.5*std::pow(sys->N1().getU()*0.1,2.0));
    deleteAll(*sys.get());
}
*/

///The following test verifies the pressure output of the ideal mixer SERJ case
/**(Figure 4.27 of <a href="http://oleg.lib.strath.ac.uk/R/?func=dbin-jump-full&object_id=26892">
 * "Toward automated design of Combined Cycle Propulsion." by Alessandro Mogavero</a>)
 */
/*TEST(CFD, PressureIdeal) {
    double tolerance = 0.01;
    std::unique_ptr<hypro::systemModule> sys = CFD(CASE1tuned,1.0);

    EXPECT_LT(1-std::fabs(sys->N1().getPress()/221273), tolerance)<<"The pressure of inlet node = " << sys->N1().getPress();
    double internalPressure[7] = {213144, 219186, 268412, 5.57961e+06, 36984.1, 298093, 166506};
    for (unsigned i = 0; i<sys->nodes_.size();i++) {
        EXPECT_LT(1-std::fabs(sys->nodes_[i].get()->getPress()/internalPressure[i]), tolerance)<<"The pressure of node " << i << " = " << sys->nodes_[i].get()->getPress();
    }
    EXPECT_LT(1-std::fabs(sys->N2().getPress()/21038.6), tolerance)<<"The pressure of exit node = " << sys->N2().getPress();

    deleteAll(*sys.get());
}
*/

///The following test verifies the pressure output of the 70% efficiency mixer SERJ case
/**(Figure 4.27 of <a href="http://oleg.lib.strath.ac.uk/R/?func=dbin-jump-full&object_id=26892">
 * "Toward automated design of Combined Cycle Propulsion." by Alessandro Mogavero</a>)
 */
/*TEST(CFD, Pressure70pct) {
    double tolerance = 0.01;
    std::unique_ptr<hypro::systemModule> sys = CFD(CASE1tuned,0.700);

    EXPECT_LT(1-std::fabs(sys->N1().getPress()/221273), tolerance)<<"The pressure of inlet node = " << sys->N1().getPress();
    double internalPressure[7] = {213144, 219186, 227630, 5.57961e+06, 36984.1, 263446, 148521};
    for (unsigned i = 0; i<sys->nodes_.size();i++) {
        EXPECT_LT(1-std::fabs(sys->nodes_[i].get()->getPress()/internalPressure[i]), tolerance)<<"The pressure of node " << i << " = " << sys->nodes_[i].get()->getPress()
                            ;
    }
    EXPECT_LT(1-std::fabs(sys->N2().getPress()/22743.8), tolerance)<<"The pressure of exit node = " << sys->N2().getPress();

    deleteAll(*sys.get());
}
*/

///The following test verifies the Mach output of the ideal mixer SERJ case
/**(Figure 4.28 of <a href="http://oleg.lib.strath.ac.uk/R/?func=dbin-jump-full&object_id=26892">
 * "Toward automated design of Combined Cycle Propulsion." by Alessandro Mogavero</a>)
 */
/*TEST(CFD, MachIdeal) {
    double tolerance = 0.01;
    std::unique_ptr<hypro::systemModule> sys = CFD(CASE1tuned,1.0);

    EXPECT_LT(std::fabs(1-sys->N1().M()/0.294509), tolerance)<<"The Mach of inlet node = " << sys->N1().M();
    double internalMach[7] = {0.377273, 0.31743, 0.461347, 1.00, 4.185, 0.225733, 1.00};
    for (unsigned i = 0; i<sys->nodes_.size();i++) {
        EXPECT_LT(std::fabs(1-sys->nodes_[i].get()->M()/internalMach[i]), tolerance)<<"The Mach of node " << i << " = " << sys->nodes_[i].get()->M();
    }
    EXPECT_LT(std::fabs(1-sys->N2().M()/2.38027), tolerance)<<"The Mach of exit node = " << sys->N2().M();

    deleteAll(*sys.get());
}
*/

///The following test verifies the Mach output of the 70% efficiency mixer SERJ case
/**(Figure 4.28 of <a href="http://oleg.lib.strath.ac.uk/R/?func=dbin-jump-full&object_id=26892">
 * "Toward automated design of Combined Cycle Propulsion." by Alessandro Mogavero</a>)
 */
/*TEST(CFD, Mach70pct) {
    double tolerance = 0.01;
    std::unique_ptr<hypro::systemModule> sys = CFD(CASE1tuned,0.700);

    EXPECT_LT(std::fabs(1-sys->N1().M()/0.294509), tolerance)<<"The Mach of inlet node = " << sys->N1().M();
    double internalMach[7] = {0.377273, 0.31743, 0.540502, 1.00, 4.185, 0.255119, 1.00};
    for (unsigned i = 0; i<sys->nodes_.size();i++) {
        EXPECT_LT(std::fabs(1-sys->nodes_[i].get()->M()/internalMach[i]), tolerance)<<"The Mach of node " << i << " = " << sys->nodes_[i].get()->M();
    }
    EXPECT_LT(std::fabs(1-sys->N2().M()/2.26279), tolerance)<<"The Mach of exit node = " << sys->N2().M();

    deleteAll(*sys.get());
}
*/
