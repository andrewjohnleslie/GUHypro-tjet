/*!
 *  \brief			Deprecated class
 *  \Deprecated
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

#ifndef TERMMODULE_H_
#define TERMMODULE_H_

#include "core/propModule.h"

namespace hypro {

	class TermModule : public propModule {
		Node ownedN_;

	public:
		TermModule();
		TermModule(Node *N2);

		TermModule(std::string name, Node *N2, int ID = defID_);

		~TermModule();

		propModuleSERIAL(TermModule);

	virtual bool canChoke()const;

		virtual GPObject &duplicate() { return *(new TermModule(*this)); }

		double calculate();

		double drag() const;
	};
}
#endif /* TERMMODULE_H_ */
