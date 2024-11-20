/*!
 *  \author     Alessandro Mogavero
 *  \n
 *  \email      alessandro.mogavero@strath.ac.uk
 *  \version    1.0
 *  \copyright  Copyright 2017 Alessandro Mogavero
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

#include "ParameterDefinitionGene.h"
#include <propModule.h>
#include <ParMixer.h>
#include <nonLinear.h>
#include <MyPopulation.h>
#include <GPTypedGeneSet.h>

namespace hypro{
ParameterDefinitionGene::ParameterDefinitionGene (double par, int value)
:
	NoNodeGene(value),
	par_(par){
}

ParameterDefinitionGene::ParameterDefinitionGene(const ParameterDefinitionGene& gpo)
:
		NoNodeGene(gpo),
		par_(gpo.par_){
}

MyGene& ParameterDefinitionGene::operator=(const MyGene& gpo){
	NoNodeGene::operator =(gpo);
	par_ = ((ParameterDefinitionGene&)gpo).par_;
	return *this;
}

ParameterDefinitionGene::~ParameterDefinitionGene(){
}

void ParameterDefinitionGene::createSystem(OptimSystem& sys, systemModule::ModulePtr parent)const{

}

double ParameterDefinitionGene::operate(double min, double max)const{
	return min + (max - min)*par_;
}

// Grow GPs according to the given parameters
void ParameterDefinitionGene::create (enum GPCreationType ctype, int allowableDepth,
		GPTypedGeneSet& ns)
{
	//Do nothing. The ParameterDefinitionGene is always a terminal
}

int ParameterDefinitionGene::isA (){
	return MyPopulation::ParameterDefinitionGeneID;
}

MyGene::Type ParameterDefinitionGene::type()const{
	return MyGene::PARAMETER;
}

std::vector<MyGene::Type> ParameterDefinitionGene::argsType()const{
	return std::vector<MyGene::Type>();
}

int ParameterDefinitionGene::isFunction()const{
	return false;
}

int ParameterDefinitionGene::isTerminal()const{
	return true;
}

void ParameterDefinitionGene::printMathStyle (ostream& os, int lastPrecedence){
	os << "Parameter-" << par_;
}

std::string ParameterDefinitionGene::name()const{
	return "Parameter-" + std::to_string(par_);
}
}
