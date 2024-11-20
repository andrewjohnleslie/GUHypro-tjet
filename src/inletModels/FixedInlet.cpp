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
#include "FixedInlet.h"
#include "ZitaInlet.h"
#include "solvers/EffIsenDuct.h"

namespace hypro {
	FixedInlet::FixedInlet(Node *N1, Node *N2)
			:
			FixedInlet(defName_, N1, N2) {
	}

	FixedInlet::FixedInlet(std::string name, Node *N1, Node *N2)
			:
			Inlet(name, N1, N2) {
		clear();

		initNodes(1, *N1);

		/// N1 is the free stream, node[0] is the intake (the first section area of the engine which meets the free stream) and N2 is the end of the inlet


		ModulePtr mod(new ZitaInlet("Forebody", N1_, nodes_[0].get()));
		add(mod);

		ModulePtr mod1(new EffIsenDuct("Inlet Duct", nodes_[0].get(), N2_));
		EffIsenDuct *mod1c = (EffIsenDuct *) mod1.get();
		mod1c->choked_ = true;
		add(mod1);
	}

	FixedInlet::FixedInlet() :
			Inlet(),
			pdropmax_(0) {
	}

	FixedInlet::FixedInlet(std::istream &is, const NodeMap &nodeMap, const ModuleMap &moduleMap)
			:
			Inlet(is, nodeMap, moduleMap) {
	}

	FixedInlet::~FixedInlet() {
		// TODO Auto-generated destructor stub
	}

	double FixedInlet::calculate() {

		double rChok = Inlet::calculate();

		N1_->A((N2_->mfr() / (N1_->getU() * N1_->rho())));

		return rChok;
	}


	void FixedInlet::reduce(const double &q) {

		std::list<ModulePtr>::iterator it = modules_.begin();
		++it;
		EffIsenDuct *itc = (EffIsenDuct *) (*it).get();
		itc->etap_ = 1 - q;

	}

	double FixedInlet::reduce() const {

		std::list<ModulePtr>::const_iterator it = modules_.begin();
		++it;
		EffIsenDuct *itc = (EffIsenDuct *) (*it).get();

		return 1 - itc->etap_;

	}

	void FixedInlet::unreduce() {

		std::list<ModulePtr>::iterator it = modules_.begin();
		++it;
		EffIsenDuct *itc = (EffIsenDuct *) (*it).get();
		itc->etap_ = 1 / pdropmax_;

	}

void FixedInlet::serialize(Archive& ar) const{
	Inlet::serialize(ar);

	ar.put("Maximum pressure drop", pdropmax_);
}

void FixedInlet::unserialize(const Archive& ar) {
	Inlet::unserialize(ar);

	pdropmax_ = ar.getDouble("Maximum pressure drop");
}
}
