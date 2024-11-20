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

#include "ModuleAdditionGene.h"
#include <ParMixer.h>
#include <MyPopulation.h>
#include <FeedbackDefinitionGene.h>
#include <GPTypedGeneSet.h>
#include <ParameterDefinitionGene.h>
#include <Combustion.h>
#include <Friction.h>
#include <InjectionPhi.h>


namespace hypro
{
ModuleAdditionGene::ModuleAdditionGene(GPNode &gpo)
	: MyGene(gpo), parameters_()
{
	if (((MyNode &)gpo).type_ != MyNode::MODULE)
	{
		throw std::logic_error(
			"Error: ModuleAdditionGene node must be of type MODULE.");
	}
}

ModuleAdditionGene::ModuleAdditionGene(const ModuleAdditionGene &gpo)
	: MyGene(gpo), parameters_()
{
	for (std::vector<std::unique_ptr<Parameter>>::const_iterator it =
			 gpo.parameters_.begin();
		 it != gpo.parameters_.end(); it++)
	{
		parameters_.push_back(std::unique_ptr<Parameter>((*it)->duplicate()));
	}
}

MyGene &ModuleAdditionGene::operator=(const MyGene &gpo)
{
	MyGene::operator=(gpo);

	ModuleAdditionGene &gpo1 = (ModuleAdditionGene &)gpo;

	parameters_.clear();
	for (std::vector<std::unique_ptr<Parameter>>::const_iterator it =
			 gpo1.parameters_.begin();
		 it != gpo1.parameters_.end(); it++)
	{
		parameters_.push_back(std::unique_ptr<Parameter>((*it)->duplicate()));
	}

	return *this;
}

ModuleAdditionGene::~ModuleAdditionGene() {}

void ModuleAdditionGene::createSystem(OptimSystem &sys,
									  systemModule::ModulePtr parent) const
{
	// add module
	systemModule::ModulePtr module(&(propModule &)node->duplicate());
	systemModule::NodePtr N2(new Node(module->N2()));
	module->N2(*N2);

	std::ostringstream name;

	// Run child genes
	for (int i = begin(MODULE); i < end(MODULE); i++)
	{
		// create modules and connect to input nodes
		NthMyChild(i)->createSystem(sys, module);
		systemModule::ModulePtr childModule1 = sys.modules_.back();
		module->inNode(i, childModule1->N2());
	}
	for (int i = begin(FEEDBACK); i < end(FEEDBACK); i++)
	{
		// Define feedback
		NthMyChild(i)->createSystem(sys, module);
	}
	uint j = 0;
	for (int i = begin(PARAMETER); i < end(PARAMETER); i++)
	{
		// Define parameters
		parameters_[j]->set(*module, (ParameterDefinitionGene &)*NthMyChild(i));
		j++;
	}

	// Link to free stream if required
	if (module->inputFreeStream())
	{
		systemModule::NodePtr N1(new Node(module->N1()));
		module->N1(*N1);
		sys.nodes_.push_back(N1);
		name.str("");
		name << N1->Name_ << "-" << sys.nodes_.size() - 1;
		N1->Name_ = name.str();
	}

	// Add output node to system
	sys.nodes_.push_back(N2);
	name.str("");
	name << N2->Name_ << "-" << sys.nodes_.size() - 1;
	N2->Name_ = name.str();

	// Add created module to system
	sys.add(module);
	name.str("");
	name << module->Name_ << "-" << sys.modules_.size() - 1;
	module->Name_ = name.str();
}

// Grow GPs according to the given parameters
void ModuleAdditionGene::create(enum GPCreationType ctype, int allowableDepth,
								GPTypedGeneSet &gs)
{
	// Loop through the whole function arguments and create new
	// offsprings (the function parameters).  Remember: There is already
	// a Gene, namely the object this routine is called for by
	// GP::create(), and it is a function.
	std::vector<Type> argsT(argsType());

	for (int n = 0; n < containerSize(); n++)
	{
		Type argType = argsT[n];

		int chooseTerm;

		if (argType == FEEDBACK || argType == PARAMETER)
		{
			// Feedback type is always a terminal
			chooseTerm = 1;
		}
		else
		{
			// Now decide whether the offspring should be a function or a
			// terminal.  If there is no more allowable depth, choose a
			// terminal in any case.
			if (ctype == GPGrow)
			{
				chooseTerm = 0;
			}
			else
			{
				// 50/50% chance of getting a function or a terminal
				chooseTerm = GPrand() % 2;
			}
			if (allowableDepth <= 1)
				chooseTerm = 1;
		}

		// Choose function or terminal
		MyGene *g;
		if (chooseTerm)
			g = &gs.chooseTerminal(argType);
		else
			g = &gs.chooseFunction(argType);

		// Create a new gene with the chosen node and place the new gene
		// into the container
		put(n, *g);

		// If the node was a function, call function recursively with
		// decreased allowable depth
		g->create(ctype, allowableDepth - 1, gs);
	}
}

int ModuleAdditionGene::isA() { return MyPopulation::ModuleAdditionGeneID; }

MyGene::Type ModuleAdditionGene::type() const { return MyGene::MODULE; }

MyNode::Type ModuleAdditionGene::nodeType() const { return MyNode::MODULE; }

std::vector<MyGene::Type> ModuleAdditionGene::argsType() const
{
	propModule *module = (propModule *)node;

	std::vector<MyGene::Type> argsType;
	if (!module->inputFreeStream())
	{
		argsType.assign(module->inNode().size(), MODULE);
	}
	if (module->canChoke())
	{
		argsType.push_back(FEEDBACK);
	}

	for (std::vector<std::unique_ptr<Parameter>>::const_iterator it =
			 parameters_.begin();
		 it != parameters_.end(); it++)
	{
		argsType.push_back(PARAMETER);
	}

	return argsType;
}

ModuleAdditionGene::Parameter::Parameter(double min, double max)
	: max_(max), min_(min) {}

ModuleAdditionGene::Parameter::Parameter(const Parameter &obj)
	: max_(obj.max_), min_(obj.min_) {}

ModuleAdditionGene::Parameter::~Parameter() {}

template <class Module>
ModuleAdditionGene::DefinedPar<Module>::DefinedPar(double min, double max,
												   double Module::*par,
												   int nodeI)
	: Parameter(min, max), par_(par), meth_(NULL), nodeI_(nodeI) {}

template <class Module>
ModuleAdditionGene::DefinedPar<Module>::DefinedPar(
	double min, double max, void (Module::*meth)(const double &), int nodeI)
	: Parameter(min, max), par_(NULL), meth_(meth), nodeI_(nodeI) {}

template <class Module>
ModuleAdditionGene::DefinedPar<Module>::DefinedPar(const DefinedPar &obj)
	: Parameter(obj), par_(obj.par_), meth_(obj.meth_), nodeI_(obj.nodeI_) {}

template <class Module>
ModuleAdditionGene::DefinedPar<Module>::~DefinedPar() {}

template <class Module>
void ModuleAdditionGene::DefinedPar<Module>::set(
	propModule &mod, const ParameterDefinitionGene &gene)
{
	if (par_)
	{
		((Module &)mod).*par_ = gene.operate(min_, max_);
	}
	else
	{
		(((Module &)mod).*meth_)(gene.operate(min_, max_));
	}
}

template <>
void ModuleAdditionGene::DefinedPar<Node>::set(
	propModule &mod, const ParameterDefinitionGene &gene)
{
	Node *n;
	if (nodeI_ > -1)
	{
		if (!mod.isCollection())
		{
			throw std::logic_error(
				"Error: Attempt to set internal node parameter of a non collection "
				"module.");
		}
		n = ((Collection &)mod).nodes_[nodeI_].get();
	}
	else
	{
		n = &mod.N2();
	}
	if (par_)
	{
		n->*par_ = gene.operate(min_, max_);
	}
	else
	{
		(n->*meth_)(gene.operate(min_, max_));
	}
}

template <class Module>
ModuleAdditionGene::Parameter *
ModuleAdditionGene::DefinedPar<Module>::duplicate() const
{
	return new DefinedPar<Module>(*this);
}

void ModuleAdditionGene::addParameter(Parameter *par)
{
	parameters_.push_back(std::unique_ptr<Parameter>(par));
	reserveSpace(containerSize() + 1);
}

template class ModuleAdditionGene::DefinedPar<Node>;
template class ModuleAdditionGene::DefinedPar<InjectionPhi<Friction>>;
template class ModuleAdditionGene::DefinedPar<Combustion<Friction>>;
}