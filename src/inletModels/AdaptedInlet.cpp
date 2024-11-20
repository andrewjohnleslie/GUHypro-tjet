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

#include "AdaptedInlet.h"
#include "NeutralLink.h"
#include "SupersonicInlet.h"
#include "chokedConv.h"

namespace hypro {
	AdaptedInlet::AdaptedInlet(std::string name, Node *N1, Node *N2, int ID)
			:
			Inlet(name, N1, N2, ID) {
		clear();

		initNodes(2, *N1_);

		ModulePtr mod2(new NeutralLink("Adaption flow", N1_, nodes_[0].get()));
		add(mod2);

		ModulePtr mod(new SupersonicInlet("Inlet Duct", nodes_[0].get(), N2_));
		SupersonicInlet *modc = (SupersonicInlet *) mod.get();
		modc->choked_ = true;
		add(mod);

		ModulePtr mod1(new chokedConv("Convergent", nodes_[0].get(), nodes_[1].get(), chokedConv::CALCULATED));
		add(mod1);
	}

	AdaptedInlet::AdaptedInlet() :
			Inlet() {
	}

	AdaptedInlet::AdaptedInlet(std::istream &is, const NodeMap &nodeMap, const ModuleMap &moduleMap)
			:
			Inlet(is, nodeMap, moduleMap) {
	}

	GPObject &AdaptedInlet::duplicate() {
		return Inlet::duplicate();
	}

	AdaptedInlet::~AdaptedInlet() {
		// TODO Auto-generated destructor stub
	}

	void AdaptedInlet::reduce(const double &A) {
		N1_->A(A);
	}

	double AdaptedInlet::reduce() const {
		return N1_->A();
	}

	bool AdaptedInlet::isreduceable() const {
		return true;
	}

	void AdaptedInlet::unreduce() {
		N1_->A(nodes_[0]->Amax_);
	}

}

//BOOST_CLASS_EXPORT_IMPLEMENT(AdaptedInlet);
