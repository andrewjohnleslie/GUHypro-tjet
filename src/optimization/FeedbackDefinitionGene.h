/*!
 *  \brief      GP gene implementation for HyPro
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

#ifndef FEEDBACKDEFINITIONGENE_H_
#define FEEDBACKDEFINITIONGENE_H_

#include <NoNodeGene.h>
#include <systemModule.h>

namespace hypro
{
class GPTypedGeneSet;

class FeedbackDefinitionGene : public NoNodeGene
{
	uint steps_;

  public:
	FeedbackDefinitionGene() : NoNodeGene(),
							   steps_(0) {}

	///Constructor based on GP node
	FeedbackDefinitionGene(uint steps, int value);

	///Copy constructor
	FeedbackDefinitionGene(const FeedbackDefinitionGene &gpo);

	virtual MyGene &operator=(const MyGene &gpo);

	virtual ~FeedbackDefinitionGene();

	FeedbackDefinitionGene &operator=(FeedbackDefinitionGene N);

	///Instantiate a new object of as copy of this object
	virtual GPObject &duplicate() { return *(new FeedbackDefinitionGene(*this)); }

	///Create the GP systemModel object starting from the GP tree
	/**
	 * This function is called by the root GP gene that in turns call this same function
	 * for its child and so on.
	 * So the GP tree is spanned recursively.
	 * @param sys reference to the system module being created.
	 * @param parent pointer to the module created by the parent Gene
	 */
	void createSystem(OptimSystem &sys, systemModule::ModulePtr parent = NULL) const;

	void create(enum GPCreationType ctype, int allowableDepth, GPTypedGeneSet &ns);

	///Returns the type ID of this class. Used for serialisation
	virtual int isA();

	///Instantiate a permanent object of type FeedbackDefinitionGene
	virtual GPObject *createObject() { return new FeedbackDefinitionGene(); }

	virtual MyGene::Type type() const;

	virtual std::vector<Type> argsType() const;

	virtual int isFunction() const;

	virtual int isTerminal() const;

	virtual void printMathStyle(ostream &os, int lastPrecedence = 0);

	virtual std::string name() const;
};
}
#endif /* FEEDBACKDEFINITIONGENE_H_ */
