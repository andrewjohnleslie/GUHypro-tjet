#include "Np3.h"
#include "US76.h"
#include <boost/math/tools/roots.hpp>
#include <chrono>

using namespace hypro;

void postprocessing(outputData &DataOut){


	DataOut.sys->printOut(std::cout);
	std::cout << std::endl << "Thrust is: " << DataOut.thrust << " N" << std::endl << std::endl;


}



int main(int argc, char *argv[]) {

    auto start = std::chrono::high_resolution_clock::now();


	double Mach = 0.8;
//	double altitude = 10668;
//	US76 atmos(altitude);
//	double p = atmos.p();
//	double T = atmos.T();
	double p = 23842;
	double T = 218.81;
	double throttle = 0.388;
//	double throttle = 1;

	int verbosity = 1;

	std::pair<bool, outputData> Run = Np3(p, T, Mach, throttle, verbosity);

	outputData DataOut = Run.second;

	postprocessing(DataOut);


	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	std::cout << "Time Elapsed: " << elapsed.count() << " seconds" << std::endl;


	std::cout << DataOut.sys->Name_ << " Completed" <<std::endl;



	return 0;
}
