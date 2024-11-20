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

#include "ComplexityPenalty.h"
#include <fstream>
#include <MyGene.h>
#include <MyPopulation.h>
#include <fenv.h>
#include <ConvNozzle.h>
#include <ModuleAdditionGene.h>
#include <GPTypedGeneSet.h>

namespace hypro{
ComplexityPenalty::ComplexityPenalty()
:
		CrossAreaConstraint(){
}

ComplexityPenalty::ComplexityPenalty(std::string name, Node* N1, Node* N2, int genes)
:
	CrossAreaConstraint(name,N1,N2,genes){
}

ComplexityPenalty::ComplexityPenalty(const ComplexityPenalty& module)
:
	CrossAreaConstraint(module){
}

ComplexityPenalty::~ComplexityPenalty() {
	// TODO Auto-generated destructor stub
}

GPObject& ComplexityPenalty::duplicate(){
	return *((GP*)(new ComplexityPenalty(*this)));
}

void ComplexityPenalty::calcFitness(){
	CrossAreaConstraint::calcFitness();

	if(stdFitness != std::numeric_limits<double>::infinity()){
		stdFitness *= 1.0 + 0.1*modules_.size();
	}
}

int ComplexityPenalty::isA (){
	return MyPopulation::ComplexityPenaltyID;
}
}