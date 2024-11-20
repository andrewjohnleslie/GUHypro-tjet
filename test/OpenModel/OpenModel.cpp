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

#define DO_QUOTE(X)       #X
#define QUOTE(X)          DO_QUOTE(X)

#include <iostream>
#include <fstream>

#include "core/systemModule.h"

int main(int argc, char *argv[]) {
	if(argc<2){
		throw std::runtime_error("Error: 1 input requested.");
	}

	thermoState::warning_ = false;

	//Read HyPro model
	std::string fileName(argv[1]);

	auto sys = std::dynamic_pointer_cast<systemModule>(propModule::loadJson(fileName));

	sys->unreduce();
	sys->calculate();

	sys->printOut(std::cout);
	std::cout << "Thrust = " << sys->thrust() << std::endl;

	std::vector<double> Xmfr = sys->propellantMfr();
	double mfr = 0.0;
	for(std::size_t i=0; i<Xmfr.size(); i++){
		mfr += Xmfr[i];
	}

	std::cout << "Isp = " << sys->thrust()/(mfr*9.81) << std::endl;

	sys->show(argc,argv);

	sys->deleteAll(sys.get());

	return 1;
}
