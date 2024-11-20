/*! \author     Alessandro Mogavero
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

#include "GP.h"
#include <gtest/gtest.h>

///This test the performances of the GP optimiser in subsonic conditions
/**<a href="http://oleg.lib.strath.ac.uk/R/?func=dbin-jump-full&object_id=26892">
 * "Toward automated design of Combined Cycle Propulsion." by Alessandro Mogavero</a>
 */
TEST(GP, subsonic) {
	propModule::chokingMach_ = 0.9;

	const double tolerance = 0.01;
	const double expFit = 0.000770562;

	std::string oldInDir = inDir;
	ostrstream strInFile;
	strInFile  << oldInDir << "/subsonic" << ends;
	inDir = strInFile.str();
	std::string seedFile = std::string(inDir) + "/seed.dat";

	setup(true, seedFile);
	MyPopulation* pop = run(0,false);
	const double fit = pop->NthGP(pop->bestOfPopulation)->getFitness();

	std::cout << fit << std::endl;
	EXPECT_LT(std::fabs(1.0 - fit/expFit), tolerance) << "The fitness = " << fit;

	inDir = oldInDir;

	delete pop;

	systemModule::resetLoopLimits_ = false;
	nonLinear<systemModule>::maxIter_ = 1000;
}

///This test the performances of the GP optimiser in supersonic conditions
/**<a href="http://oleg.lib.strath.ac.uk/R/?func=dbin-jump-full&object_id=26892">
 * "Toward automated design of Combined Cycle Propulsion." by Alessandro Mogavero</a>
 */
TEST(GP, supersonic) {
	propModule::chokingMach_ = 0.9;

	const double tolerance = 0.01;
	const double expFit = 0.000296763;

	std::string oldInDir = inDir;
	ostrstream strInFile;
	strInFile  << oldInDir << "/supersonic" << ends;
	inDir = strInFile.str();
	std::string seedFile = std::string(inDir) + "/seed.dat";

	setup(true, seedFile);
	MyPopulation* pop = run(0,false);
	const double fit = pop->NthGP(pop->bestOfPopulation)->getFitness();

	std::cout << fit << std::endl;
	EXPECT_LT(std::fabs(1.0 - fit/expFit), tolerance) << "The fitness = " << fit;

	inDir = oldInDir;
	delete pop;

	systemModule::resetLoopLimits_ = false;
	nonLinear<systemModule>::maxIter_ = 1000;
}
