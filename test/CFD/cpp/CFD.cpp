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

#include <boost/property_tree/json_parser.hpp>
#include "Archive.h"
#include "CFD.h"

int main(int argc, char *argv[]) {

    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
    // Case inputs
	caseTag selCase = CASE1tuned;
	double eta = 1.0;
	for(int i=1; i<argc; i+=2){
		if(strcmp(argv[i],"-case")==0){
			sscanf(argv[i+1],"%d",(int*)&selCase);
		}else if(strcmp(argv[i],"-eta")==0){
			sscanf(argv[i+1],"%f",(float*)&eta);
		}else{
			std::cerr << "Error: wrong input." << std::endl;
			exit(1);
		}
	}

	std::unique_ptr<hypro::systemModule> sys = CFD(selCase,eta);

    double k = 1.5*std::pow(sys->N1().getU()*0.1,2.0); //turbulent kinetic energy with 10% of turbulence intensity
    double omega = std::sqrt(k)/1; //turbulent specific dissipation rate with length scale=1m

	std::cout << "Inlet conditions:" << std::endl
	<< "\tgamma = " << sys->N1().gamma() << std::endl
	//<< "\tT0 = " << sys->N1().gas().TH(sys->N1().H(),600.0) << std::endl
	<< "\tflowRate = " << sys->N1().mfr()*2/360 << std::endl
	<< "\tk = " << k << std::endl
	<< "\tomega = " << omega << std::endl;
	std::cout << "Parameters:" << std::endl << "\teta = " << eta << std::endl;

	std::vector<double> Mfr = sys->propellantMfr();
	hypro::Collection::OutType out(sys->out());

    sys->printOut(std::cout);

    std::stringstream ss;
    sys->saveJson(ss);
	std::cout << ss.str() << std::endl;

	std::shared_ptr<hypro::propModule> sys1 = hypro::propModule::loadJson(ss);

	std::stringstream ss1;
	sys1->saveJson(ss1);
	std::cout << ss1.str() << std::endl;

    // create and open a character archive for output
	std::ofstream ofs("filename");
	ofs << ss1.str();
	ofs.close();

	deleteAll(*sys.get());
	return 0;
}
