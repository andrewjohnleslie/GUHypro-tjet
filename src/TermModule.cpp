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

#include "TermModule.h"

namespace hypro {
	TermModule::TermModule(Node *N2)
			:
			TermModule(defName_, N2) {
	}

	TermModule::TermModule(std::string name, Node *N2, int ID)
			:
			propModule(name, (Node *) NULL, &ownedN_, ID, 0),
			ownedN_(*N2) {
	}

	TermModule::~TermModule() {
		// TODO Auto-generated destructor stub
	}

	double TermModule::calculate() {
		if (verbosity_ > 0) {
			std::cout << "Warning: class TermModule does not need to be calculated." << std::endl;
		}
		return 1;
	}

double TermModule::drag()const{
	return N2_->A()*N2_->rho()*std::pow(N2_->getU(),2.0);
}

bool TermModule::canChoke()const{
	return false;
}
}