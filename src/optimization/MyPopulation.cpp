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

#include "MyPopulation.h"
#include <ModuleAdditionGene.h>
#include <FeedbackDefinitionGene.h>
#include <ParameterDefinitionGene.h>
#include <CrossAreaConstraint.h>
#include <ComplexityPenalty.h>

namespace hypro {
MyPopulation::MyPopulation(GPVariables& GPVar_, GPAdfNodeSet& adfNs_, Node& freeStream, Node& exit)
:
	GPPopulation (GPVar_, adfNs_),
	freeStream_(&freeStream),
	exit_(&exit),
	GPtype_(new OptimSystem()){
}

MyPopulation::MyPopulation (MyPopulation& gpo)
:
	GPPopulation(gpo),
	freeStream_(gpo.freeStream_),
	exit_(gpo.exit_),
	GPtype_((OptimSystem*)(GP*)&gpo.GPtype_->duplicate()){
}

// Prints out the complete population
	void MyPopulation::printOn(ostream &os) {
		for (int n = 0; n < containerSize(); n++) {
			OptimSystem *current = (OptimSystem *) NthGP(n);
			if (current)
				os << n << " " << *(GP *) current << "\t" << current->getFitness() << " " << current->message_ << endl;
			else
				os << "(NULL)" << endl;
		}
	}

void MyPopulation::registerOptimClasses ()
{
	GPRegisterClass (new MyPopulation());
	GPRegisterClass ((GP*)(new OptimSystem()));
	GPRegisterClass ((GP*)(new CrossAreaConstraint()));
	GPRegisterClass ((GP*)(new ComplexityPenalty()));
	GPRegisterClass (new ModuleAdditionGene());
	GPRegisterClass (new FeedbackDefinitionGene());
	GPRegisterClass (new ParameterDefinitionGene());
}

	char *MyPopulation::load(istream &is) {
		char *out = GPPopulation::load(is);

		for (int n = 0; n < containerSize(); n++) {
			OptimSystem *current = (OptimSystem *) NthGP(n);
			if (current) {
				current->N1(*freeStream_);
				current->N2(*exit_);
			}
		}
		return out;
	}

GP* MyPopulation::createGP (int numOfGenes){
	OptimSystem* out = (OptimSystem*)(GP*)GPtype_->createObject();
	out->reserveSpace(numOfGenes);
	out->N1(*freeStream_);
	out->N2(*exit_);

	return out;
}
}
