/*
// Created by Mark De Luca on 14/12/18.
//
// For more information, contact m.de-luca.1@research.gla.ac.uk
*/

#include "US76.h"
#include "TurboFan.h"
#include <boost/math/tools/roots.hpp>
#include <chrono>


using namespace hypro;

enum DP { OFF, ON };

double kH_HP {0};
double kH_LP {0};
double T04_ {0};
double T045_ {0};
double Aratio_HP {0};
double Aratio_LP {0};


void run_model (const int verbosity, const vector<double> AltitudeRange, const vector<double> MachRange,
        const vector<double> throttleRange, const double AltitudeDP, const double MachDP, const double throttleDP, plotData ReturnedOutputs[]){

		std::size_t counter {0};
		int ondesign;
		plotData output;
		offDesignData Run;
		double CaptureArea {0};


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



				for (std::size_t ithrottle = 0; ithrottle<(throttleRange.size()); ithrottle++) {


//					std::cout << "Starting New Throttle Position! " << "\n" << "----------------------------------------" << "\n\n";

					const double throttle = throttleRange[ithrottle];

					std::cout << "\n\n********************************\n\n\n";
					fmt::print(std::cout, "{:<12} {:<12.2f}\n","Mach Number", Mach);
					fmt::print(std::cout, "{:<12} {:<12.2f}\n", "Altitude ", Altitude);
					fmt::print(std::cout, "{:<12} {:<12.4f}\n\n", "Throttle", throttle);


					if (throttle == throttleDP && Altitude == AltitudeDP && Mach == MachDP)
					{

						std::cout << "Is On Design" << "\n";
						ondesign = ON;

						std::pair<bool, outputData> DP = TurboFan (pDP, TDP, MachDP, throttleDP, 1, ondesign, CaptureArea, T04_, T045_, Aratio_HP, Aratio_LP, kH_HP, kH_LP);

						Run.analysis(output, DP.second, DP.second);
						ReturnedOutputs[counter] = output;
//		                DP.second.sys->saveJson("/home/mark/Documents/MATLAB/OffDesign/pyCycleComparison/ConcordeSST/DUMPpycycleHPLPHyProRFPDP.txt");


					}
					else
					{

						CaptureArea = 0;
						std::cout << "Is Off Design" << "\n";
						std::cout << "First Running at Design Point" << "\n\n";
						ondesign = ON;


						std::pair<bool, outputData> DP = TurboFan (pDP, TDP, MachDP, throttleDP, verbosity, ondesign, CaptureArea, T04_, T045_, Aratio_HP, Aratio_LP, kH_HP, kH_LP);
//		                DP.second.sys->saveJson("/home/mark/Documents/MATLAB/OffDesign/pycycleHPLPHyProDP.txt");


						std::cout << "\n" << "- - - - - - - - - - - -  - - - - - - - - - - - - - - - - - -";
						std::cout << "\n\n" << "Now Running Off Design" << "\n\n";
						ondesign = OFF;

						T04_ = DP.second.T04;
						kH_HP = DP.second.kH_HP;
						Aratio_HP = DP.second.Aratio_HP;
						T045_ = DP.second.T045;
						kH_LP = DP.second.kH_LP;
						Aratio_LP = DP.second.Aratio_LP;


						Run.analysis(output, DP.second, DP.second);

						if (counter == 0){
							ReturnedOutputs[counter] = output;
							counter++;
						}

						std::pair<bool, outputData> offDP;
						offDP = TurboFan (p, T, Mach, throttle, verbosity, ondesign, CaptureArea, T04_, T045_, Aratio_HP, Aratio_LP, kH_HP, kH_LP);
//

						Run.analysis(output, DP.second, offDP.second);
						CaptureArea = (output.minod/offDP.second.min);
//						ReturnedOutputs[counter] = output;

						std::cout <<  "Rerunning Off Design with Corrected Capture Area...\n\n" ;
						offDP = TurboFan (p, T, Mach, throttle, verbosity, ondesign, CaptureArea, T04_, T045_, Aratio_HP, Aratio_LP, kH_HP, kH_LP);


//						Run.analysis(output, DP.second, offDP.second);
//						output.FT = offDP.second.thrust;
						ReturnedOutputs[counter] = output;
//		                offDP.second.sys->saveJson("/home/mark/Documents/MATLAB/OffDesign/pycycleHPLPHyProODT0.txt");


					}

					counter++;

				}

			}

		}

}




int main(int argc, char *argv[]) {

	auto start = std::chrono::high_resolution_clock::now();

	const int verbosity = 1;
	const int cppaste = 0;


	std::vector<double> MachRange = {0.3241};
	std::vector<double> AltitudeRange = {5000};
	std::vector<double> throttleRange = {0.4};

	const double MachDP = 0.5;
	const double AltitudeDP = 5000;
	const double throttleDP = 0.4;

//	std::vector<double> MachRange = {MachDP};
//	std::vector<double> AltitudeRange = {AltitudeDP};
//	std::vector<double> throttleRange = {throttleDP};

	int numbOfPoints;
	if (((MachRange.size() *  AltitudeRange.size() * throttleRange.size()) == 1) && (throttleRange.front() == throttleDP && AltitudeRange.front() == AltitudeDP && MachRange.front() == MachDP)){
		numbOfPoints = 1;
	}
	else{
		numbOfPoints = 1 + (MachRange.size() *  AltitudeRange.size() * throttleRange.size());
	}

	plotData ReturnedOutputs[numbOfPoints];
	run_model(verbosity, AltitudeRange, MachRange, throttleRange, AltitudeDP, MachDP, throttleDP, ReturnedOutputs);

	offDesignData Run;
	Run.disp(ReturnedOutputs, numbOfPoints, cppaste);

	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	std::cout << "\n" << "Time Elapsed: " << elapsed.count() << " seconds" << "\n";


	return 0;


}
