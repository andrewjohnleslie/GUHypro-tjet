/*!
// Created by Mark De Luca on 07/12/18.
//
// For any queries, contact m.de-luca.1@research.gla.ac.uk
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
 * along with HyPro.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "TurboJet.h"
#include "US76.h"
#include <boost/math/tools/roots.hpp>
#include <chrono>
#include <iostream>
#include <time.h>

using namespace hypro;

enum DP { OFF, ON };

double kH_ {0};
double T04_ {0};
double Aratio_ {0};



void run_model (const int verbosity, const vector<double> AltitudeRange, const vector<double> MachRange,
        const vector<double> throttleRange, const double AltitudeDP, const double MachDP, const double throttleDP, plotData ReturnedOutputs[]){


		std::size_t counter = 0;
		int ondesign;
		plotData output;
		offDesignData Run;
		double CaptureArea {0};
		
		//std::cout << "TESTING - TURBOJET.CPP - LINE 49 " << std::endl; 	// TRIGGERED


		for (std::size_t iMach = 0; iMach<(MachRange.size()); iMach++) {

//			std::cout << "Starting New Mach Range! " << "\n" << "----------------------------------------" << "\n\n" ;
			const double Mach = MachRange[iMach];
			if (Mach <= 0.0){
				throw std::out_of_range("Please select a non zero (positive) number for Mach");
			}


			for (std::size_t iAlt = 0; iAlt<(AltitudeRange.size()); iAlt++) {

//				std::cout << "Starting New Altitude Range! " << "\n" << "----------------------------------------" << "\n\n";

				const double Altitude = AltitudeRange[iAlt];

				US76 atmo = US76(Altitude);
				US76 atmoDP = US76(AltitudeDP);

				double p = atmo.p();
				double T = atmo.T();
				double pDP = atmoDP.p();
				double TDP = atmoDP.T();

				//std::cout << "TESTING - TURBOJET.CPP - LINE 75 " << std::endl;	// TRIGGERED
				
				for (std::size_t ithrottle = 0; ithrottle<(throttleRange.size()); ithrottle++) {

//					std::cout << "New Throttle Position! " << "\n" << "----------------------------------------" << "\n\n";

					const double throttle = throttleRange[ithrottle];


					fmt::print(std::cout, "{:<12} {:<12.2f}\n","\nMach Number", Mach);
					fmt::print(std::cout, "{:<12} {:<12.2f}\n", "Altitude ", Altitude);
					fmt::print(std::cout, "{:<12} {:<12.4f}\n\n", "Throttle", throttle);


					if (throttle == throttleDP && Altitude == AltitudeDP && Mach == MachDP)
					{
						std::cout << "Is On Design" << "\n";
						ondesign = ON;

						std::cout << "... Calling Engine Function \n";
						std::cout << "_ _ _ _ _ _ _ _ _ _\n\n\n";
						//std::cout << "TESTING - TURBOJET.CPP - LINE 97 " << std::endl;	// TRIGGERED
						std::pair<bool, outputData> DP     = TurboJet (p, T, Mach, throttle, verbosity, ondesign, CaptureArea, T04_, Aratio_, kH_);

						//std::cout << "TESTING - TURBOJET.CPP - LINE 100 " << std::endl;	// TRIGGERED
						
						Run.analysis(output, DP.second, DP.second);
						ReturnedOutputs[counter] = output;
		                DP.second.sys->saveJson("/home/andrew/Desktop/JsonTurboJet.txt");

					}
					else
					{
						//std::cout << "TESTING - TURBOJET.CPP - LINE 109 " << std::endl;	// TRIGGERED
						CaptureArea = 0;
						std::cout << "Is Off Design" << "\n";
						std::cout << "First Running at Design Point..." << "\n\n";
						ondesign = ON;

						std::pair<bool, outputData> DP     = TurboJet (pDP, TDP, MachDP, throttleDP, verbosity, ondesign, CaptureArea, T04_, Aratio_, kH_);

						std::cout << "\n" << "- - - - - - - - - - - -  - - - - - - - - - - - - - - - - - -";
						std::cout << "\n\n" << "Now Running Off Design..." << "\n\n";
						ondesign = OFF;
						//std::cout << "TESTING - TURBOJET.CPP - LINE 120 " << std::endl;	// TRIGGERED
						T04_ = DP.second.T04;
						kH_ = DP.second.kH;
						Aratio_ = DP.second.Aratio;

						Run.analysis(output, DP.second, DP.second);
						//std::cout << "TESTING - TURBOJET.CPP - LINE 126 " << std::endl;	// TRIGGERED
						if (counter == 0){
							ReturnedOutputs[counter] = output;
							counter++;
						}

						//std::pair<bool, outputData> offDP;
						//std::cout << "TESTING - TURBOJET.CPP - LINE 133 " << std::endl;	// TRIGGERED
						std::pair<bool, outputData> offDP  = TurboJet (p,   T,   Mach,   throttle,   verbosity, ondesign, CaptureArea, T04_, Aratio_, kH_);
						//std::cout << "TESTING - TURBOJET.CPP - LINE 135 " << std::endl;	// NOT TRIGGERED

						Run.analysis(output, DP.second, offDP.second);
						CaptureArea = (output.minod/offDP.second.min);


						std::cout <<  "Rerunning Off Design with Corrected Capture Area...\n\n" ;
						//std::cout << "TESTING - TURBOJET.CPP - LINE 142 " << std::endl;	// 
//						
						offDP  = TurboJet (p,   T,   Mach,   throttle,   verbosity, ondesign, CaptureArea, T04_, Aratio_, kH_);
//						Run.analysis(output, DP.second, offDP.second);	// Should this be un-commented
//						output.FT = offDP.second.thrust;			// Should this be un-commented
						ReturnedOutputs[counter] = output;
		                offDP.second.sys->saveJson("/home/andrew/Desktop/JsonTurboJet.txt");


					}

					counter++;

				}

			}

		}

}




int main(int argc, char *argv[]) {
			
		//std::cout << "TESTING - TURBOJET.CPP - LINE 168 " << std::endl;	// TRIGGERED
		auto start = std::chrono::high_resolution_clock::now();

		const int verbosity = 0;	// Gives greater detail of testcase [SET TO 1]
		const int cppaste = 1;		// Provides the post-processing value print (throttle, mfr_in, T04/T02, pi_c, m2A2, mbar) [SET TO 1]
		std::vector<double> AltitudeRange = {0};		// TESTING
		std::vector<double> MachRange = {0.5};		// TESTING
		//std::vector<double> MachRange = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2};		// TESTING
		std::vector<double> throttleRange = {0.25995};	// TESTING
//		std::vector<double> AltitudeRange = {0};
//		std::vector<double> AltitudeRange = {1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000, 13000, 14000, 15000, 16000, 17000, 18000};
		//std::vector<double> AltitudeRange = {2000, 4000, 6000, 8000, 10000, 11000, 12000, 14000, 16000, 18000};
//		std::vector<double> MachRange = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0};
//		std::vector<double> MachRange = {0.01};
//		std::vector<double> throttleRange = {0.194};
//		std::vector<double> throttleRange = {0.21775};
//		std::vector<double> throttleRange = {0.21774, 0.217741,	0.217742, 0.217743,	0.217744, 0.217745,	0.217746,
//				0.217747,	0.217748, 0.217749, 0.217751, 0.217752, 0.217753, 0.217754, 0.217755, 0.217756,0.217757, 0.217758, 0.217759, 0.21776};
//		std::vector<double> throttleRange = {0.21, 0.211, 0.212, 0.213, 0.214, 0.215, 0.216, 0.217001, 0.218, 0.21901, 0.22};
//		std::vector<double> throttleRange = {0.17, 0.175, 0.18, 0.185, 0.19, 0.195, 0.2, 0.205, 0.21, 0.215, 0.22, 0.225, 0.23, 0.235, 0.24};
		int numbOfPoints = 1 + (MachRange.size() *  AltitudeRange.size() * throttleRange.size());

		const double MachDP		= 0.01;
		const double AltitudeDP = 0;
		const double throttleDP = 0.25995; //0.21775;
//		const double throttleDP = 0.194;
		//std::cout << "TESTING - TURBOJET.CPP - LINE 189 " << std::endl;	// TRIGGERED
		plotData ReturnedOutputs[numbOfPoints];
		//std::cout << "TESTING - TURBOJET.CPP - LINE 191 " << std::endl;	// TRIGGERED
		run_model(verbosity, AltitudeRange, MachRange, throttleRange, AltitudeDP, MachDP, throttleDP, ReturnedOutputs);
		//std::cout << "TESTING - TURBOJET.CPP - LINE 193 " << std::endl;	// TRIGGERED
		offDesignData Run;
		Run.disp(ReturnedOutputs, numbOfPoints, cppaste);

		auto finish = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = finish - start;
		std::cout << "\nTime Elapsed: " << elapsed.count() << " seconds" << "\n";

		return 0;

}
