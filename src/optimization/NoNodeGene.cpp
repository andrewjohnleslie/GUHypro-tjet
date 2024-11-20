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

#include "NoNodeGene.h"
#include <propModule.h>
#include <ParMixer.h>
#include <nonLinear.h>
#include <MyPopulation.h>
#include <GPTypedGeneSet.h>

namespace hypro {
NoNodeGene::NoNodeGene (int value)
:
	MyGene(),
	value_(value){
}

NoNodeGene::NoNodeGene(const NoNodeGene& gpo)
:
		MyGene(gpo),
		value_(gpo.value_){
}

MyGene& NoNodeGene::operator= (const MyGene& gpo){
	MyGene::operator=(gpo);
	value_ = ((NoNodeGene&)gpo).value_;
	return *this;
}

NoNodeGene::~NoNodeGene(){
}

MyNode::Type NoNodeGene::nodeType()const{
	throw std::logic_error("Error: This Gene type has not a Node.");
}

propModule& NoNodeGene::getNode()const{
	throw std::logic_error("Error: This Gene type has not a Node.");
}

int NoNodeGene::value()const{
	return value_;
}

std::string NoNodeGene::name()const{
	throw std::logic_error("Error: Method name needs to be implemented.");
}

char* NoNodeGene::load (istream& is)
{
  is >> value_;

  // Load container
  return GPContainer::load (is);
}
}