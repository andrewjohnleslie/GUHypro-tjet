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

#include "AdaptedThroatInlet.h"
#include "NeutralLink.h"
#include "chokedConv.h"
#include "SupersonicInlet.h"
#include "ModuleGraph.h"

namespace hypro {
	AdaptedThroatInlet::AdaptedThroatInlet(std::string name, Node *N1, Node *N2, int ID)
			:
			Inlet(name, N1, N2, ID) {
		clear();

		initNodes(2, *N1_);

		ModulePtr mod1(new chokedConv("Convergent", N1_, nodes_[1].get(), chokedConv::ASSIGNED));
		add(mod1);

		ModulePtr mod2(new NeutralLink("Adaption flow", N1_, nodes_[0].get()));
		add(mod2);

		ModulePtr mod(new SupersonicInlet("Inlet Duct", nodes_[0].get(), N2_));
		add(mod);
	}

	AdaptedThroatInlet::AdaptedThroatInlet()
			:
			Inlet() {
	}

	GPObject &AdaptedThroatInlet::duplicate() {
		AdaptedThroatInlet *newMod(new AdaptedThroatInlet(*this));
		duplicate(newMod);

		return *newMod;
	}

	AdaptedThroatInlet::AdaptedThroatInlet(std::istream &is, const NodeMap &nodeMap, const ModuleMap &moduleMap)
			:
			Inlet(is, nodeMap, moduleMap) {
	}

	double AdaptedThroatInlet::calculate() {
		double rChok;
		try {
			rChok = Inlet::calculate();
		} catch (std::range_error &e) {
			if (N1_->M() > 1.0) {
				N1_->A(nodes_[0]->Amax_);
				chokedConv *mod1c = (chokedConv *) modules_.front().get();
				mod1c->throat_ = chokedConv::CALCULATED;

				rChok = Inlet::calculate();

				mod1c->throat_ = chokedConv::ASSIGNED;
			} else {
				nodes_[0]->A(nodes_[0]->Amax_);

				list<ModulePtr>::iterator it = modules_.begin();
				ModulePtr mod2 = *(++it);
				ModulePtr mod2new(new isentropicDuct("Adaption flow", N1_, nodes_[0].get()));
				exchange(mod2, mod2new);

				rChok = Collection::calculate();

				exchange(mod2new, mod2);
			}
		}

		return rChok;
	}

	AdaptedThroatInlet::~AdaptedThroatInlet() {
		// TODO Auto-generated destructor stub
	}

	void AdaptedThroatInlet::reduce(const double &A) {
		nodes_[1]->A(A);
	}

	double AdaptedThroatInlet::reduce() const {
		return nodes_[1]->A();
	}

	bool AdaptedThroatInlet::isreduceable() const {
		return true;
	}

	void AdaptedThroatInlet::unreduce() {
		nodes_[1]->A(nodes_[1]->Amax_);
	}

	Glib::RefPtr<ModuleGraph> AdaptedThroatInlet::draw() {
		graph_ = ModuleGraph::create(this);
		return graph_;
	}
}

//BOOST_CLASS_EXPORT_IMPLEMENT(AdaptedThroatInlet);
