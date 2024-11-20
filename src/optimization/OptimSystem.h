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

#ifndef OPTIMSYSTEM_H_
#define OPTIMSYSTEM_H_

#include "core/systemModule.h"
#include <core/FeedBack.h>
#include <examplewindow.h>
#include <gp.h>
#include <gtkmm.h>

namespace hypro
{
class MyGene;
class GPTypedGeneSet;

class OptimSystem : public systemModule, public GP
{

	static const int swapAttempts_ = 5; ///<Maximum number of attempt performed during mutation operation

  protected:
	///Calculates the fitness once the model has been run
	/**Called by the fitness() method after the HyPro model has
	 * been fully calculated.
	 * This can be redefined in subclasses in order to imlement different fitness functions.
	 */
	virtual void calcFitness();

  public:
	static int printTexStyle;   ///<1 if the GP genome should be print in tex style. 0 otherwise
	static bool checkValidity_; ///<True if the validity of the genome should be checked at runtime

	std::string message_; ///<Message recorded at the end of fitness evaluation

	///Constructor with no inputs
	OptimSystem();

	///Constructor
	/**
	 * @param name system module name
	 * @param N1 pointer to system input node
	 * @param N2 pointer to system output node
	 * @param genes total number of genes
	 */
	OptimSystem(std::string name, Node *N1, Node *N2, int genes);

	///Copy constructor
	OptimSystem(const OptimSystem &module);

	~OptimSystem();

	///Duplicate object
	virtual GPObject &duplicate();

	///Create a gene
	/**
	 * Used in class GP to create the root gene of the GP.
	 * It creates a ModuleAdditionGene because the root gene cannot be otherwise.
	 *
	 * @param gpo Gene's underlying GP node
	 * @return POinter to the new Gene
	 */
	virtual GPGene *createGene(GPNode &gpo);

	///Print the object
	virtual void printOn(ostream &os);

	///Print the object on standard output
	void printOn();

	///Returns a pointer to the nth gene
	MyGene *NthMyGene(int n);

	///Returns a pointer to the nth gene
	MyGene *NthMyGene(int n) const;

	///Evaluate the fitness function of GP system.
	/**It first create the systemModule calling the createSystem method of the root gene.
	 * Then it calculate the systemModule and eventually the fitness function
	 * calling method calcFitness().
	 */
	virtual void evaluate();

	/// Crossover operator
	/** Cross the objects contained in the given container.  The function
	 * is responsible for the given container and has to delete it, if a
	 * completely new container is returned.  We use only the objects in
	 * the container for the crossover operation, and not the object for
	 * which we were called (though this is in our case a member of the
	 * container).  We can use all container objects and do with them what
	 * we want.  We only have to return another container.  The returned
	 * container needn't contain any objects, but definitely should at
	 * least sometimes, otherwise an infinite loop will occur, because the
	 * generate-function tries to fill the new population with new
	 * members, and if we return none, well...
	 * @param parents GP breeding parents
	 * @param maxdepthforcrossover maximum depth of the generated offspring
	 * @return pointer to the container containing the offspring
	 */
	GPContainer &cross(GPContainer *parents, int maxdepthforcrossover);

	///Mutate the GP with the probability given by parameter in the GPVariables.
	/**
	 * @param GPVar GP setting
	 * @param adfNs node set to be used to create mutation
	 */
	void mutate(GPVariables &GPVar, GPAdfNodeSet &adfNs);

	///Shrink mutation operator
	/**This block of code performs shrink mutation on a genetic program.
	 * A function node is chosen by random from a random GP tree, and one
	 * of the children of the function takes the position of the parent.
	 */
	void shrinkMutation();

	///Swap mutation operator
	/**This block of code perform allele or swap mutation on a genetic
	 * program.  This comprised of code being swapped with other code with
	 * certain constraints.  Any terminal can be swapped with any other
	 * but functions can only be swapped with other functions with the
	 * same arguments.  This means that the mutation does not have to
	 * create new branches when different function types are swapped which
	 * seems implicitly wrong
	 * @param adfNs node set to be used to create mutation
	 */
	void swapMutation(GPAdfNodeSet &adfNs);

	void create(enum GPCreationType ctype, int allowableDepth,
				GPAdfNodeSet &adfNs);

	///Check validity of the GP genome
	void checkValidity() const;

	///Return class ID. Used in serialisation.
	virtual int isA();

	///Instantiate new OptimSystem object
	/**createObject does not initialise the N1_ and N2_, they have to be assigned manually after loading!
	 * @return pointer to created object
	 */
	virtual GPObject *createObject() { return (GP *)(new OptimSystem()); }

	///Load GP genome from stream
	virtual char *load(istream &is);

	///Creates GUI tree representing the GP genome
	int createTree(int argc, char *argv[]);

	void resolveNodeValues(GPAdfNodeSet &adfNs);
};
}
#endif /* OPTIMSYSTEM_H_ */
