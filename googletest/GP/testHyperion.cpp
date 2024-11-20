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
//testHyperion.cpp -- tests Hyperion example case (thrust and fuel mass flow)
#include "US76.h"
#include "Hyperion.h"
#include <gtest/gtest.h>
//#include "fmt/format.h"

///This test is to ensure that the thrust predicted by the Hyperion Ejector case agrees with
/**Figure 4.9 of <a href="http://oleg.lib.strath.ac.uk/R/?func=dbin-jump-full&object_id=26892">
 * "Toward automated design of Combined Cycle Propulsion." by Alessandro Mogavero</a>
 */
TEST(HyperionEjector, Thrust) {
    double tolerance = 0.01;

    int Machiterations = 7;
    int Altiterations = 4;

    double MachRange[7] = {0.01, 0.5, 1.0, 1.5, 2.0, 2.5, 2.9};
    double AltitudeRange[4] = {0, 6096, 12192, 18288};

    double expectedThrust[4][9] = {{92.3134,94.1007,104.157,118.941,186.173,311.509,529.395},
                                   {85.3985,87.388,98.1248,109.64,140.048,210.644,315.6},
                                   {80.5147,81.8224,87.951,96.1152,117.364,147.455,189.136},
                                   {79.5474,80.152,82.8729,85.9673,96.0348,112.211,131.076}};

    for(int j=0; j<Altiterations; j++){
        double Altitude_ = AltitudeRange[j];
        US76 atmo(Altitude_);
        for(int i=0; i<Machiterations; i++){
            double Mach = MachRange[i];
            double Thrust;
            std::vector<double> Mfr;
            Collection::OutType out;
            // Hyperion(atmo.p(),atmo.T(),Mach,1,0,0,Thrust,Mfr,out, dynamic_pressure, emissions, total_flow);
            systemModule* sys = Hyperion(atmo.p(),atmo.T(),Mach,0,0,Thrust,Mfr,out);
            fmt::print("Reference: {} vs Actual: {}\n", expectedThrust[j][i], Thrust/4.449717900000000e+03);
            EXPECT_LT(1-((Thrust/(4.449717900000000e+03))/expectedThrust[j][i]), tolerance);
            systemModule::deleteAll(sys);
        }
    }
}

 ///This test is to ensure that the ISP predicted by the Hyperion Ejector case agrees with
 /**Figure 4.10 of <a href="http://oleg.lib.strath.ac.uk/R/?func=dbin-jump-full&object_id=26892">
  * "Toward automated design of Combined Cycle Propulsion." by Alessandro Mogavero</a>
  */
 TEST(HyperionEjector, ISP) {
     double tolerance = 0.01;

     int Machiterations = 7;
     int Altiterations = 4;

     double MachRange[7] = {0.01, 0.5, 1.0, 1.5, 2.0, 2.5, 2.9};
     double AltitudeRange[4] = {0, 6096, 12192, 18288};

     double expectedISP[4][7] = {{408.976,414.085,451.226,495.679,708.069,1038.36,1413.71},
                                 {386.753,394.41,437.455,477.301,586.713,809.099,1075.21},
                                 {369.207,374.643,400.535,432.509,516.031,625.812,763.946},
                                 {366.928,369.504,381.252,393.657,435.783,500.682,572.039}};

     for (int j = 0; j < Altiterations; j++) {
         double Altitude_ = AltitudeRange[j];
         US76 atmo(Altitude_);
         for (int i = 0; i < Machiterations; i++) {
             double Mach = MachRange[i];
             double Thrust;
             std::vector<double> Mfr;
             Collection::OutType out;
             systemModule* sys = Hyperion(atmo.p(),atmo.T(),Mach,0,0,Thrust,Mfr,out);
             // Hyperion(atmo.p(),atmo.T(),Mach,1,0,0,Thrust,Mfr,out, dynamic_pressure, emissions, total_flow);
             EXPECT_LT(1-((Thrust/(9.81*(Mfr[0]+Mfr[1]+Mfr[2]+Mfr[3])))/expectedISP[j][i]), tolerance);
             systemModule::deleteAll(sys);
         }
     }
 }

 ///This test is to ensure that the Coefficient of thrust predicted by the pure Hyperion Ramjet case agrees with
 /**Figure 4.14 of <a href="http://oleg.lib.strath.ac.uk/R/?func=dbin-jump-full&object_id=26892">
  * "Toward automated design of Combined Cycle Propulsion." by Alessandro Mogavero</a>
  */
 TEST(HyperionPureRamjet, ThrustCoefficient) {
     double tolerance = 0.01;

     int Machiterations = 9;
     int Altiterations = 4;

     double MachRange[9] = {2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0};
     double AltitudeRange[4] = {9144, 15240, 21336, 33528};


     double expectedCT[4][9] = {{0.284366,0.506062,0.806913,1.18168,1.51571,1.36564,1.21331,1.07529,0.953257},
                                 {0.284162,0.504978,0.803642,1.17549,1.57766,1.42578,1.26876,1.1257,0.999165},
                                 {0.284183,0.505095,0.803997,1.17616,1.57055,1.41887,1.26239,1.11991,0.993861},
                                 {0.284444,0.506463,0.808093,1.18393,1.49518,1.34583,1.19511,1.05867,0.938156}};

     for (int j = 0; j < Altiterations; j++) {
         double Altitude_ = AltitudeRange[j];
         US76 atmo(Altitude_);
         for (int i = 0; i < Machiterations; i++) {
             double Mach = MachRange[i];
             double Thrust;
             std::vector<double> Mfr;
             Collection::OutType out;
             systemModule *sys = Hyperion(atmo.p(),atmo.T(),Mach,1,0,Thrust,Mfr,out);
             // Hyperion(atmo.p(),atmo.T(),Mach,1,0,Thrust,Mfr,out, dynamic_pressure, emissions, total_flow);
             EXPECT_LT(1-((Thrust/(0.5*sys->N1().rho()*std::pow(sys->N1().getU(), 2)*2.5084))/expectedCT[j][i]), tolerance);
             systemModule::deleteAll(sys);
         }
     }
 }

 ///This test is to ensure that the ISP predicted by the pure Hyperion Ramjet case agrees with
 /**Figure 4.15 of <a href="http://oleg.lib.strath.ac.uk/R/?func=dbin-jump-full&object_id=26892">
  * "Toward automated design of Combined Cycle Propulsion." by Alessandro Mogavero</a>
  */
 TEST(HyperionPureRamjet, ISP) {
     double tolerance = 0.01;

     int Machiterations = 9;
     int Altiterations = 4;

     double MachRange[9] = {2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0};
     double AltitudeRange[4] = {9144, 15240, 21336, 33528};


     double expectedISP[4][9] = {{2493.2, 3201.3, 3407.38, 3257.56, 3081.91, 3137.31, 3107.18, 3036.65, 2942.6},
                                 {2576.38, 3308.8, 3521.92, 3370.52, 3122.13, 3187.92, 3162.46, 3094.57, 3002.57},
                                 {2566.85, 3296.49, 3508.74, 3357.49, 3117.64, 3182.23, 3156.25, 3088.1, 2995.79},
                                 {2465.71, 3165.75, 3369.59, 3220.23, 3068.04, 3120.15, 3088.43, 3016.96, 2922.33}};

     for (int j = 0; j < Altiterations; j++) {
         double Altitude_ = AltitudeRange[j];
         US76 atmo(Altitude_);
         for (int i = 0; i < Machiterations; i++) {
             double Mach = MachRange[i];
             double Thrust;
             std::vector<double> Mfr;
             Collection::OutType out;
             systemModule* sys = Hyperion(atmo.p(),atmo.T(),Mach,1,0,Thrust,Mfr,out);
             // Hyperion(atmo.p(),atmo.T(),Mach,1,0,Thrust,Mfr,out, dynamic_pressure, emissions, total_flow);
             EXPECT_LT(1-((Thrust / (9.81 * Mfr[0]))/expectedISP[j][i]), tolerance);
             systemModule::deleteAll(sys);
         }
     }
 }
