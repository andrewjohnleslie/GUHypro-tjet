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
#include "../pyCycleTjetComparison.h"

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

//				double p = atmo.p();
//				double T = atmo.T();
				double p = 84307.056852;
				double T = 278.243888888889;
//				double pDP = atmoDP.p();
//				double TDP = atmoDP.T();
//				double p = 101324.703484;
//				double T = 288.15;
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

						std::pair<bool, outputData> DP     = pyCycleTjetComparison (pDP, TDP, Mach, throttle, verbosity, ondesign, T04_, Aratio_, kH_, CaptureArea);


						Run.analysis(output, DP.second, DP.second);

						ReturnedOutputs[counter] = output;
//		                DP.second.sys->saveJson("/home/mark/Documents/MATLAB/OffDesign/HyProTwoMethodCoincidenceCheckDP.txt");



					}



					else
					{

						std::cout << "Is Off Design" << std::endl;
						std::cout << "First Running at Design Point" << std::endl << std::endl;
						ondesign = ON;
						CaptureArea = 0;

						std::pair<bool, outputData> DP     = pyCycleTjetComparison (pDP, TDP, MachDP, throttleDP, verbosity, ondesign, T04_, Aratio_, kH_, CaptureArea);

						T04_ = DP.second.T04;
						kH_ = DP.second.kH;
						Aratio_ = DP.second.Aratio;
						std::cout  << std::endl << "T04_ is " << T04_ << std::endl;
						std::cout  << std::endl << "kH_ is " << kH_ << std::endl;
						std::cout << std::endl << "- - - - - - - - - - - -  - - - - - - - - - - - - - - - - - -";
						std::cout << std::endl << std::endl << "Now Running Off Design" << std::endl;
						ondesign = OFF;
						CaptureArea = MachDP/Mach;

						if (counter == 0){

							Run.analysis(output, DP.second, DP.second);
							ReturnedOutputs[counter] = output;
							counter++;

						}
//		                DP.second.sys->saveJson("/home/mark/Documents/MATLAB/OffDesign/pycycleTjetHyProDPComparisonMarch2021.txt");


						std::pair<bool, outputData> offDP;
						offDP = pyCycleTjetComparison (p,   T,   Mach,   throttle,  verbosity, ondesign, T04_, Aratio_, kH_, CaptureArea);
						kH_ = offDP.second.kH;
						std::cout  << std::endl << "kH_ is " << kH_ << std::endl;


						Run.analysis(output, DP.second, offDP.second);
						CaptureArea = CaptureArea * (output.minod/offDP.second.min);

//						std::cout <<  "Rerunning Off Design with Corrected Capture Area...\n\n" ;
//						offDP  = pyCycleTjetComparison (p,   T,   Mach, throttle, verbosity	, ondesign, T04_, Aratio_, kH_, CaptureArea);
//
//						Run.analysis(output, DP.second, offDP.second);
//						CaptureArea = CaptureArea * (output.minod/offDP.second.min);
//
//
//						std::cout <<  "Rerunning Off Design with Corrected Capture Area...\n\n" ;
//						offDP  = pyCycleTjetComparison (p,   T,   Mach, throttle, verbosity	, ondesign, T04_, Aratio_, kH_, CaptureArea);
//
//
//						Run.analysis(output, DP.second, offDP.second);
						ReturnedOutputs[counter] = output;

//						ostrstream SaveString;
//						SaveString << "/home/mark/Documents/MATLAB/OffDesign/OffDesignAnalytical/HyProJsonDump/TjetOffDP-throttle-" << throttle << "-Mach-" << Mach << "-Alt-" << Altitude << ".txt" << std::ends;
//						std::ofstream fil(SaveString.str());
//		                offDP.second.sys->saveJson(fil);


		                offDP.second.sys->saveJson("/home/mark/Documents/MATLAB/OffDesign/pycycleTjetHyProOD2ComparisonMarch2021.txt");


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
		std::vector<double> AltitudeRange = {1524.5};
		std::vector<double> MachRange = {1.2};
		std::vector<double> throttleRange = {0.1985};

//		std::vector<double> AltitudeRange = {0};
//		std::vector<double> MachRange = {0.000001};
//		std::vector<double> throttleRange = {0.25995};
//		std::vector<double> throttleRange = {0.246635};


//		std::vector<double> throttleRange = {0.18, 0.185, 0.19, 0.195, 0.2, 0.205, 0.21, 0.215, 0.22, 0.225, 0.23, 0.235, 0.24, 0.245, 0.25, 0.255, 0.26, 0.265, 0.27, 0.275, 0.28};
//		std::vector<double> throttleRange = {0.255, 0.2555, 0.256, 0.2565, 0.257, 0.2575, 0.258, 0.2585, 0.259, 0.2595, 0.26, 0.2605, 0.261, 0.2615, 0.262, 0.2625, 0.263, 0.2635, 0.264, 0.2645, 0.265};
//		std::vector<double> throttleRange = {0.2599, 0.25991, 0.25992, 0.25993, 0.25994, 0.25996, 0.25997, 0.25998, 0.25999};
//		std::vector<double> MachRange = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6 ,0.7, 0.8, 0.9, 1.001, 1.1, 1.2, 1.3, 1.4, 1.6, 1.8, 2.0};
//		std::vector<double> MachRange = {1.6, 1.8, 2.0};
//		std::vector<double> AltitudeRange = {1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000, 13000, 14000, 15000, 16000, 17000, 18000};
//		std::vector<double> AltitudeRange = {8000};



		int numbOfPoints = 1 + (MachRange.size() *  AltitudeRange.size() * throttleRange.size());

		const double MachDP = 0.000001;
//		const double MachDP = 0.01;
		const double AltitudeDP = 0;
		const double throttleDP = 0.25995;



		plotData ReturnedOutputs[numbOfPoints];
		run_model(verbosity, AltitudeRange, MachRange, throttleRange, AltitudeDP, MachDP, throttleDP, ReturnedOutputs);

		offDesignData Run;
		Run.disp(ReturnedOutputs, numbOfPoints, cppaste);

		auto finish = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = finish - start;
		std::cout << "Time Elapsed: " << elapsed.count() << " seconds" << std::endl;

		return 0;

}
