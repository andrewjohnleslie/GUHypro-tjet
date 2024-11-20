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

#ifndef MODULEADDITIONGENE_H_
#define MODULEADDITIONGENE_H_

#include <MyGene.h>
#include <systemModule.h>

namespace hypro
{
class ParameterDefinitionGene;

class ModuleAdditionGene : public MyGene
{
	class Parameter
	{
	  protected:
		double max_;
		double min_;

	  public:
		Parameter(double min, double max);
		Parameter(const Parameter &obj);
		virtual ~Parameter();

		virtual Parameter *duplicate() const = 0;

		virtual void set(propModule &mod, const ParameterDefinitionGene &gene) = 0;
	};

	std::vector<std::unique_ptr<Parameter>> parameters_;

  public:
	template <class Module>
	class DefinedPar : public Parameter
	{
		double Module::*par_;
		void (Module::*meth_)(const double &);
		int nodeI_;

	  public:
		DefinedPar(double min, double max, double Module::*par, int nodeI = -1);
		DefinedPar(double min, double max, void (Module::*meth)(const double &),
				   int nodeI = -1);
		DefinedPar(const DefinedPar &obj);
		virtual ~DefinedPar();
		virtual Parameter *duplicate() const;

		void set(propModule &mod, const ParameterDefinitionGene &gene);
	};

	ModuleAdditionGene() : MyGene() {}

	/// Constructor based on GP node
	ModuleAdditionGene(GPNode &gpo);

	/// Copy constructor
	ModuleAdditionGene(const ModuleAdditionGene &gpo);

	virtual MyGene &operator=(const MyGene &gpo);

	virtual ~ModuleAdditionGene();

	ModuleAdditionGene &operator=(ModuleAdditionGene N);

	/// Instantiate a new object of as copy of this object
	virtual GPObject &duplicate() { return *(new ModuleAdditionGene(*this)); }

	/// Create the GP systemModel object starting from the GP tree
	/**
   * This function is called by the root GP gene that in turns call this same
   * function
   * for its child and so on.
   * So the GP tree is spanned recursively.
   * @param sys reference to the system module being created.
   * @param parent pointer to the module created by the parent Gene
   */
	void createSystem(OptimSystem &sys,
					  systemModule::ModulePtr parent = NULL) const;

	void create(enum GPCreationType ctype, int allowableDepth,
				GPTypedGeneSet &ns);

	/// Returns the type ID of this class. Used for serialisation
	virtual int isA();

	/// Instantiate a permanent object of type ModuleAdditionGene
	virtual GPObject *createObject() { return new ModuleAdditionGene(); }

	virtual MyGene::Type type() const;

	virtual MyNode::Type nodeType() const;

	virtual std::vector<Type> argsType() const;

	void addParameter(Parameter *par);
};
}
#endif /* MODULEADDITIONGENE_H_ */
