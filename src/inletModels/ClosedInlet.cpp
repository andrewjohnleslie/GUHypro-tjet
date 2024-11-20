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

#include "ClosedInlet.h"

namespace hypro {
	ClosedInlet::ClosedInlet(Node *N1, Node *N2)
			:
			propModule(defName_, N1, N2) {
	}

	ClosedInlet::ClosedInlet(std::string name, Node *N1, Node *N2)
			:
			propModule(name, N1, N2) {
	}

	ClosedInlet::ClosedInlet(std::istream &is, const NodeMap &nodeMap, const ModuleMap &moduleMap)
			:
			propModule(is, nodeMap, moduleMap) {
	}

	ClosedInlet::~ClosedInlet() {
		// TODO Auto-generated destructor stub
	}

	double ClosedInlet::calculate() {
		N1_->A(0.0);

	return 1;
}

bool ClosedInlet::canChoke()const{
	return false;
}
}
