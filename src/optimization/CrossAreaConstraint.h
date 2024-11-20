/*!
 *  \brief      Implementation of a constraint on the total cross section of the engine
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

#ifndef CROSSAREACONSTRAINT_H_
#define CROSSAREACONSTRAINT_H_

#include <OptimSystem.h>

namespace hypro{
class CrossAreaConstraint: public OptimSystem {

protected:
	///Calculates the fitness once the model has been run
	/**Called by the fitness() method after the HyPro model has
	 * been fully calculated.
	 * This can be redefined in subclasses in order to imlement different fitness functions.
	 */
	virtual void calcFitness();

public:

	///Constructor with no inputs
	CrossAreaConstraint();

	///Constructor
	/**
	 * @param name system module name
	 * @param N1 pointer to system input node
	 * @param N2 pointer to system output node
	 * @param genes total number of genes
	 */
	CrossAreaConstraint(std::string name, Node* N1, Node* N2, int genes);

	///Copy constructor
	CrossAreaConstraint(const CrossAreaConstraint& module);

	~CrossAreaConstraint();

	///Duplicate object
	virtual GPObject& duplicate ();

	///Return class ID. Used in serialisation.
	virtual int isA ();

	///Instantiate new CrossAreaConstraint object
	/**createObject does not initialise the N1_ and N2_, they have to be assigned manually after loading!
	 * @return pointer to created object
	 */
	virtual GPObject* createObject() { return (GP*)(new CrossAreaConstraint()); }
};
}
#endif /* CROSSAREACONSTRAINT_H_ */
