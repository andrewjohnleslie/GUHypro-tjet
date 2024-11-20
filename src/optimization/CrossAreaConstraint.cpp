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

#include "CrossAreaConstraint.h"
#include <fstream>
#include <MyGene.h>
#include <MyPopulation.h>
#include <fenv.h>
#include <ConvNozzle.h>
#include <ModuleAdditionGene.h>
#include <GPTypedGeneSet.h>

namespace hypro{
CrossAreaConstraint::CrossAreaConstraint()
:
		OptimSystem(){
}

CrossAreaConstraint::CrossAreaConstraint(std::string name, Node* N1, Node* N2, int genes)
:
	OptimSystem(name,N1,N2,genes){
}

CrossAreaConstraint::CrossAreaConstraint(const CrossAreaConstraint& module)
:
	OptimSystem(module){
}

CrossAreaConstraint::~CrossAreaConstraint() {
	// TODO Auto-generated destructor stub
}

GPObject& CrossAreaConstraint::duplicate(){
	return *((GP*)(new CrossAreaConstraint(*this)));
}

void CrossAreaConstraint::calcFitness(){
	OptimSystem::calcFitness();

	double S = crossArea();

	if(stdFitness != std::numeric_limits<double>::infinity()){
		double Ct_fit = 0.0;
		const double S_max = 6.0;
		if(S > S_max){
			Ct_fit = 100.0*S;
		}

		stdFitness += Ct_fit;
	}
}

int CrossAreaConstraint::isA (){
	return MyPopulation::CrossAreaConstraintID;
}
}