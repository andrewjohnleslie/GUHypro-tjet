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

#include "GPTypedGeneSet.h"
#include <MyPopulation.h>

namespace hypro{
GPTypedGeneSet::GPTypedGeneSet ()
:
	GPNodeSet(),
	genes_(),
	funcByType_(),
	termByType_(){
}

GPTypedGeneSet::GPTypedGeneSet(int n):
	GPNodeSet(n),
	genes_(),
	funcByType_(),
	termByType_(){
}

GPTypedGeneSet::GPTypedGeneSet (const GPTypedGeneSet& gpo){
	throw std::logic_error("Error: copy not allowed for class GPTypedGeneSet");
}

int GPTypedGeneSet::isA(){
	throw std::logic_error("Error: method not implemented.");
//	return MyPopulation::GPTypedGeneSetID;
}

// Put a node in the container.
void GPTypedGeneSet::putGene (MyGene* gpo)
{
	genes_.push_back(std::unique_ptr<MyGene>(gpo));

	//MyNode& gpoT = (MyNode&)gpo;

	if(gpo->isTerminal()){
		if(termByType_.find(gpo->type()) == termByType_.end()){
			termByType_.emplace(gpo->type(), std::vector<MyGene*>());
		}
		termByType_.at(gpo->type()).push_back(gpo);
	}else{
		if(funcByType_.find(gpo->type()) == funcByType_.end()){
			funcByType_.emplace(gpo->type(), std::vector<MyGene*>());
		}
		funcByType_.at(gpo->type()).push_back(gpo);
	}
}

MyGene& GPTypedGeneSet::chooseFunction(MyGene::Type type){
	std::vector<MyGene*> genesType = funcByType_.at(type);
	std::size_t numGenes = genesType.size();

	return (MyGene&)genesType.at(GPrand() % numGenes)->duplicate();
}

MyGene& GPTypedGeneSet::chooseTerminal(MyGene::Type type){
	std::vector<MyGene*> genesType = termByType_.at(type);
	std::size_t numGenes = genesType.size();

	return (MyGene&)genesType.at(GPrand() % numGenes)->duplicate();
}

MyGene* GPTypedGeneSet::chooseGeneWithArgs (int args, MyGene::Type type)
{
	std::size_t i, num;
	MyGene *n;

	std::vector<MyGene*> genesType;
	if(args == 0){
		genesType = termByType_.at(type);
	}else{
		genesType = funcByType_.at(type);
	}

	// First we count all nodes that have the specified number of
	// arguments
	for (i=0, num=0; i<genesType.size(); i++)
	if ((n=genesType.at(i)))
		if (n->containerSize()==args)
	num++;

	// No node with given number of arguments?
	if (num==0)
		return NULL;

	// Choose one
	std::size_t k=GPrand() % num;

	// Return the node with chosen index
	for (i=0, num=0; i<genesType.size(); i++)
	if ((n=genesType.at(i)))
		if (n->containerSize()==args)
	if (num++==k)
		return (MyGene*)&n->duplicate();

#if GPINTERNALCHECK
	GPExitSystem ("GPNodeSet::chooseNodeWithArgs", "Internal error");
#endif

	// Avoid compiler warnings (this code is never reached)
	return NULL;
}

MyGene* GPTypedGeneSet::searchForGene (int value){
	// Move through the container
	for (std::list<std::unique_ptr<MyGene>>::iterator it=genes_.begin();
		  it!=genes_.end(); it++){
		if ((*it)->value() == value){
			return it->get();
		}
	}

  // Not found
  return NULL;
}

}