/*!
 *  \brief      Injection model
 *  \details    Adds the injected fuel simply increasing the pressure
 *              in order to account for the increased mass flow.
 *              To be used only for replicating particular CFD setup (e.g. Lorrain test case).
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

#ifndef INJECTION_H_
#define INJECTION_H_

#include "core/propModule.h"
#include <iostream>

namespace hypro {
	class injection :
			public propModule {


public:
	std::vector<double> Xf_; ///<Composition of the injected fuel
	
	injection();

		injection(Node *N1, Node *N2);


		injection(std::string name, Node *N1, Node *N2);

		injection(const injection &mod);

		~injection();

		propModuleSERIAL(injection);

		double calculate();

	virtual bool canChoke()const;

	virtual void serialize(Archive& ar) const;
	virtual void unserialize(const Archive& ar);
};
}

#endif /* INJECTION_H_ */
