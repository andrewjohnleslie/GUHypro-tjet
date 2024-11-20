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

#include "NeutralLink.h"

namespace hypro {
	NeutralLink::NeutralLink(std::string name, Node *N1, Node *N2, int ID)
			:
			propModule(name, N1, N2, ID) {
	}

	NeutralLink::NeutralLink() :
			propModule() {
	}

	NeutralLink::NeutralLink(std::istream &is, const NodeMap &nodeMap, const ModuleMap &moduleMap)
			:
			propModule(is, nodeMap, moduleMap) {
	}

	NeutralLink::NeutralLink(const NeutralLink &mod)
			:
			propModule(mod) {
	}

	NeutralLink::~NeutralLink() {
		// TODO Auto-generated destructor stub
	}

	double NeutralLink::calculate() {
		N2_->calculate_stagnation = false;
		//copy everything with exception of Amax_ and Amin_
		//N2_->p_ = N1_->p_;
		//N2_->T_ = N1_->T_;
		//N2_->U_ = N1_->U_;
		N2_->setTP(N1_->getTemp(), N1_->getPress());
		N2_->setU(N1_->getU());
		N2_->A(N1_->A());
		N2_->X(N1_->X());

		return 1.0; //never chocked
	}


bool NeutralLink::canChoke()const{
	return false;
}
}

//BOOST_CLASS_EXPORT_IMPLEMENT(NeutralLink);
