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

#include "US76.h"
#include <boost/math/tools/roots.hpp>
#include <chrono>
#include "../pyCycleTjetMatching.h"

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


		for (std::size_t iMach = 0; iMach<(MachRange.size()); iMach++) {

			std::cout << "Starting New Mach Range! " << std::endl << "----------------------------------------" << std::endl << std::endl ;
			const double Mach = MachRange[iMach];
			if (Mach <= 0.0){
				throw std::out_of_range("Please select a non zero (positive) number for Mach");
			}


			for (std::size_t iAlt = 0; iAlt<(AltitudeRange.size()); iAlt++) {

				std::cout << "Starting New Altitude Range! " << std::endl << "----------------------------------------" << std::endl << std::endl;

				const double Altitude = AltitudeRange[iAlt];

				US76 atmo = US76(Altitude);
				US76 atmoDP = US76(AltitudeDP);

//				double p = 101324.703484;
//				double T = 288.15;
				double p = 84307.056852;
				double T = 278.243888888889;
//				double pDP = atmoDP.p();
//				double TDP = atmoDP.T();
				double pDP = 101324.703484;
				double TDP = 288.15;



				for (std::size_t ithrottle = 0; ithrottle<(throttleRange.size()); ithrottle++) {

					std::cout << "New Throttle Position! " << std::endl << "----------------------------------------" << std::endl << std::endl;

					const double throttle = throttleRange[ithrottle];



					std::cout << "Mach Number is: " << Mach << std::endl;
					std::cout << "Altitude (in metres) is: " << Altitude << std::endl;
					std::cout << "Throttle position is: " << throttle << std::endl << std::endl;


					if (throttle == throttleDP && Altitude == AltitudeDP && Mach == MachDP)
					{

						std::cout << "Is On Design" << std::endl;
						ondesign = ON;

						std::pair<bool, outputData> DP     = pyCycleTjetMatching (pDP, TDP, Mach, throttle, verbosity, ondesign, T04_, Aratio_, kH_, CaptureArea);


						Run.analysis(output, DP.second, DP.second);

						ReturnedOutputs[counter] = output;
//		                DP.second.sys->saveJson("/home/mark/Documents/MATLAB/OffDesign/pycycleTjetHyProDPMatchedMarch2021.txt");



					}



					else
					{

						std::cout << "Is Off Design" << std::endl;
						std::cout << "First Running at Design Point" << std::endl << std::endl;
						ondesign = ON;
						CaptureArea = 0;

						std::cout << "TESTING - pycycletjetmatch.cpp - LINE 119 " << std::endl;
						std::pair<bool, outputData> DP     = pyCycleTjetMatching (pDP, TDP, MachDP, throttleDP, verbosity, ondesign, T04_, Aratio_, kH_, CaptureArea);
						std::cout << "TESTING - pycycletjetmatch.cpp - LINE 121 " << std::endl;

						T04_ = DP.second.T04;
						kH_ = DP.second.kH;
						Aratio_ = DP.second.Aratio;
						std::cout  << std::endl << "T04_ is " << T04_ << std::endl;
						std::cout  << std::endl << "kH_ is " << kH_ << std::endl;
						std::cout << std::endl << "- - - - - - - - - - - -  - - - - - - - - - - - - - - - - - -";
						std::cout << std::endl << std::endl << "Now Running Off Design" << std::endl;
						ondesign = OFF;

						if (counter == 0){

							Run.analysis(output, DP.second, DP.second);
							ReturnedOutputs[counter] = output;
							counter++;

						}
//		                DP.second.sys->saveJson("/home/mark/Documents/MATLAB/OffDesign/pycycleTjetHyProDPMatchedMarch2021.txt");


						std::pair<bool, outputData> offDP;
						offDP = pyCycleTjetMatching (p,   T,   Mach,   throttle,   1, ondesign, T04_, Aratio_, kH_, CaptureArea);
						kH_ = offDP.second.kH;

						std::cout  << std::endl << "kH_ is " << kH_ << std::endl;

						Run.analysis(output, DP.second, offDP.second);
//						CaptureArea = (output.minod/offDP.second.min);
//						ReturnedOutputs[counter] = output;

//						std::cout <<  "Rerunning Off Design with Corrected Capture Area...\n\n" ;
//						offDP  = pyCycleTjetMatching (p,   T,   Mach,   throttle,   verbosity, ondesign, T04_, Aratio_, kH_, CaptureArea);
//						Run.analysis(output, DP.second, offDP.second);
//						output.minod = offDP.second.min;
						ReturnedOutputs[counter] = output;
//		                offDP.second.sys->saveJson("/home/mark/Documents/MATLAB/OffDesign/pycycleTjetHyProOD2Matched.txt");


					}

					counter++;

				}

			}

		}

}




int main(int argc, char *argv[]) {


		auto start = std::chrono::high_resolution_clock::now();

		const int verbosity = 0;
		const int cppaste = 0;
		std::vector<double> AltitudeRange = {1524};
		std::vector<double> MachRange = {1.2};
		std::vector<double> throttleRange = {0.1985};
//		std::vector<double> AltitudeRange = {0};
//		std::vector<double> MachRange = {0.000001};
//		std::vector<double> throttleRange = {0.259962466};
//		std::vector<double> throttleRange = {0.246635};
		int numbOfPoints = 1 + (MachRange.size() *  AltitudeRange.size() * throttleRange.size());

		const double MachDP = 0.000001;
		const double AltitudeDP = 0;
		const double throttleDP = 0.259962466;


		plotData ReturnedOutputs[numbOfPoints];
		run_model(verbosity, AltitudeRange, MachRange, throttleRange, AltitudeDP, MachDP, throttleDP, ReturnedOutputs);

		offDesignData Run;
		Run.disp(ReturnedOutputs, numbOfPoints, cppaste);

		auto finish = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = finish - start;
		std::cout << "Time Elapsed: " << elapsed.count() << " seconds" << std::endl;

		return 0;

}
