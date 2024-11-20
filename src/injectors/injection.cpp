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

#include "Archive.h"
#include "injection.h"

namespace hypro {
	injection::injection(Node *N1, Node *N2)
			:
			propModule(defName_, N1, N2) {
	}

	injection::injection(std::string name, Node *N1, Node *N2)
			:
			propModule(name, N1, N2) {
	}

	injection::injection(const injection &mod)
			:
			propModule(mod) {
	}

	injection::~injection() {
		// TODO Auto-generated destructor stub
	}

	double injection::calculate() {
//		std::cout << N1_->report() << std::endl;
		N2_->X(N1_->X() * (1 - Xf_[0]) + Xf_);

	//N2_->p_ = N1_->p_/(1 - Xf_[0]);
	//N2_->U_ = N1_->U_;
	//N2_->T_ = N1_->T_;
	N2_->setTP(N1_->getTemp(),N1_->getPress()/(1-Xf_[0]));
	N2_->setU(N1_->getU());
	return 1.0; //never choked
}

bool injection::canChoke()const{
	return false;
}

void injection::serialize(Archive& ar) const{
	propModule::serialize(ar);

	ar.put("Composition of injected fuel", Xf_);
}

void injection::unserialize(const Archive& ar) {
	propModule::unserialize(ar);

	Xf_ = ar.getVector("Composition of injected fuel");
}
}

