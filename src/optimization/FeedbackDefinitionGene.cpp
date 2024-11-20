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

#include "FeedbackDefinitionGene.h"
#include <propModule.h>
#include <ParMixer.h>
#include <nonLinear.h>
#include <MyPopulation.h>
#include <GPTypedGeneSet.h>


namespace hypro
{
FeedbackDefinitionGene::FeedbackDefinitionGene(uint steps, int value)
	: NoNodeGene(value),
	  steps_(steps)
{
}

FeedbackDefinitionGene::FeedbackDefinitionGene(const FeedbackDefinitionGene &gpo)
	: NoNodeGene(gpo),
	  steps_(gpo.steps_)
{
}

MyGene &FeedbackDefinitionGene::operator=(const MyGene &gpo)
{
	MyGene::operator=(gpo);
	steps_ = ((FeedbackDefinitionGene &)gpo).steps_;
	return *this;
}

FeedbackDefinitionGene::~FeedbackDefinitionGene()
{
}

void FeedbackDefinitionGene::createSystem(OptimSystem &sys, systemModule::ModulePtr parent) const
{

	//Create feedback link
	if (!sys.modules_.empty())
	{
		uint i = 1;
		if (steps_ == 0)
		{
			if (parent->isreduceable())
			{
				parent->chokingFeedback_ = parent;
			}
			else
			{
				for (Collection::TreeIterator it = --sys.end(); !it.isPastBegin(); it--)
				{
					if (it->isreduceable())
					{
						parent->chokingFeedback_ = it.getModulePtr();
						break;
					}
				}
			}
		}
		else
		{
			systemModule::ModulePtr feed;
			bool found = false;
			for (Collection::TreeIterator it = --sys.end(); !it.isPastBegin(); it--)
			{
				if (it->isreduceable())
				{
					feed = it.getModulePtr();
					if (i == steps_)
					{
						parent->chokingFeedback_ = feed;
						found = true;
						break;
					}
					i++;
				}
			}
			if (!found)
			{
				parent->chokingFeedback_ = feed;
			}
		}
	}
}

// Grow GPs according to the given parameters
void FeedbackDefinitionGene::create(enum GPCreationType ctype, int allowableDepth,
									GPTypedGeneSet &ns)
{
	//Do nothing. The FeedbackDefinitionGene is always a terminal
}

int FeedbackDefinitionGene::isA()
{
	return MyPopulation::FeedbackDefinitionGeneID;
}

MyGene::Type FeedbackDefinitionGene::type() const
{
	return MyGene::FEEDBACK;
}

std::vector<MyGene::Type> FeedbackDefinitionGene::argsType() const
{
	return std::vector<MyGene::Type>();
}

int FeedbackDefinitionGene::isFunction() const
{
	return false;
}

int FeedbackDefinitionGene::isTerminal() const
{
	return true;
}

void FeedbackDefinitionGene::printMathStyle(ostream &os, int lastPrecedence)
{
	os << "Feedback-" << steps_;
}

std::string FeedbackDefinitionGene::name() const
{
	return "Feedback-" + std::to_string(steps_);
}
}