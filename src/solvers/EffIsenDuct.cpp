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
#include "EffIsenDuct.h"

namespace hypro {
	EffIsenDuct::EffIsenDuct(std::string name, Node *N1, Node *N2, int ID)
			:
			isentropicDuct(name, N1, N2, ID),
			etap_(1.0),
			etaT_(1.0) {
	}

	EffIsenDuct::EffIsenDuct(std::istream &is, const NodeMap &nodeMap, const ModuleMap &moduleMap)
			:
			isentropicDuct(is, nodeMap, moduleMap) {
		is >> etap_ >> etaT_;
	}

	EffIsenDuct::EffIsenDuct()
			:
			isentropicDuct(),
			etap_(1.0),
			etaT_(1.0) {
	}

	EffIsenDuct::EffIsenDuct(const EffIsenDuct &mod)
			:
			isentropicDuct(mod),
			etap_(mod.etap_),
			etaT_(mod.etaT_) {
	}

	EffIsenDuct::~EffIsenDuct() {
		// TODO Auto-generated destructor stub
	}

	void EffIsenDuct::pressLoss() {
		p02_p01_ = etap_;
	}

	void EffIsenDuct::thermLoss() {
		T02_T01_ = etaT_;
	}

void EffIsenDuct::serialize(Archive& ar) const{
	isentropicDuct::serialize(ar);

	ar.put("Tot pressure ratio", etap_);
	ar.put("Tot temperature ratio", etaT_);
}

void EffIsenDuct::unserialize(const Archive& ar) {
	isentropicDuct::unserialize(ar);

	etap_ = ar.getDouble("Tot pressure ratio");
	etaT_ = ar.getDouble("Tot temperature ratio");
}
}

