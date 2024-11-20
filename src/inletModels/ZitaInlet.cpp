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
#include "ZitaInlet.h"
#include "interpolation.h"

namespace hypro {
	ZitaInlet::ZitaInlet(Node *N1, Node *N2)
			:
			ZitaInlet(defName_, N1, N2) {
	}

	ZitaInlet::ZitaInlet(std::string name, Node *N1, Node *N2)
			:
			EffMachIsenDuct(name, N1, N2) {
	}

	ZitaInlet::ZitaInlet()
			:
			EffMachIsenDuct(),
			zita_() {
	}

	ZitaInlet::ZitaInlet(std::istream &is, const NodeMap &nodeMap, const ModuleMap &moduleMap)
			:
			EffMachIsenDuct(is, nodeMap, moduleMap) {
	}

	ZitaInlet::~ZitaInlet() {
		// TODO Auto-generated destructor stub
	}

	double ZitaInlet::calculate() {
		double zita = interpolation::linear(zita_[0], zita_[1], N1_->M());
		N1_->A(zita * N1_->Amax_);
		return EffMachIsenDuct::calculate();
	}


void ZitaInlet::serialize(Archive& ar) const{
	EffMachIsenDuct::serialize(ar);

	ar.put("Ratio of critical flow", zita_);
}

void ZitaInlet::unserialize(const Archive& ar) {
	EffMachIsenDuct::unserialize(ar);

	zita_ = ar.getLookup("Ratio of critical flow");
}
}
