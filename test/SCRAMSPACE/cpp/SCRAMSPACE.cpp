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

#include "SCRAMSPACE.h"

int main(int argc, char *argv[]) {
	double p = 4100;
	double T = 370;
	double U = 2830;

	double Thrust;
	std::vector<double> Mfr;
	Collection::OutType out;

//    SCRAMSPACE(p,T,U,0.4,1,Thrust,Mfr,out);
    hypro::systemModule* sys = SCRAMSPACE(p,T,U,0.4,1,Thrust,Mfr,out);
    hypro::systemModule::deleteAll(sys);
	std::cout  << std::endl <<"Thrust = " << Thrust << ", Mfr = " << Mfr[0] <<  std::endl;

	return 0;
}
