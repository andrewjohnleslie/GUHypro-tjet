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
#include "EffMachIsenDuct.h"
//#include <boost/serialization/vector.hpp>

namespace hypro {
	EffMachIsenDuct::EffMachIsenDuct(Node *N1, Node *N2)
			:
			EffMachIsenDuct(defName_, N1, N2) {
	}

	EffMachIsenDuct::EffMachIsenDuct(std::string name, Node *N1, Node *N2)
			:
			isentropicDuct(name, N1, N2) {
	}

	EffMachIsenDuct::EffMachIsenDuct(std::istream &is, const NodeMap &nodeMap, const ModuleMap &moduleMap)
			:
			isentropicDuct(is, nodeMap, moduleMap) {
	}

	EffMachIsenDuct::EffMachIsenDuct()
			:
			isentropicDuct(),
			etap_(),
			etaT_() {
	}

	EffMachIsenDuct::EffMachIsenDuct(const EffMachIsenDuct &mod)
			:
			isentropicDuct(mod),
			etap_(mod.etap_),
			etaT_(mod.etaT_) {
	}

	EffMachIsenDuct::~EffMachIsenDuct() {
		// TODO Auto-generated destructor stub
	}

	void EffMachIsenDuct::pressLoss() {
		if (etap_[0].size()) {
			p02_p01_ = interpolation::linear(etap_[0], etap_[1], N1_->M());
		} else {
			p02_p01_ = 1.0;
		}
	}

	void EffMachIsenDuct::thermLoss() {
		if (etaT_[0].size()) {
			T02_T01_ = interpolation::linear(etaT_[0], etaT_[1], N1_->M());
		} else {
			T02_T01_ = 1.0;
		}
	}

void EffMachIsenDuct::serialize(Archive& ar) const{
	isentropicDuct::serialize(ar);

	ar.put("Tot Pressure ratio vs Mach", etap_);
	ar.put("Tot Temperature ratio vs Mach", etaT_);
}

void EffMachIsenDuct::unserialize(const Archive& ar) {
	isentropicDuct::unserialize(ar);

	etap_ = ar.getLookup("Tot Pressure ratio vs Mach");
	etaT_ = ar.getLookup("Tot Temperature ratio vs Mach");
}

}
