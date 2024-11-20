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

#include "Inlet.h"
#include "solvers/isentropicDuct.h"
#include "chokedConv.h"
#include "ModuleGraph.h"

namespace hypro {
	Inlet::Inlet(std::string name, Node *N1, Node *N2, int ID)
			:
			Collection(name, N1, N2, ID, 1) {
		initNodes(2, *N1_);

		ModulePtr mod2(new isentropicDuct("Adaption flow", N1_, nodes_[0].get()));
		add(mod2);

		ModulePtr mod(new isentropicDuct("Inlet Duct", nodes_[0].get(), N2_));
		isentropicDuct *modc = (isentropicDuct *) mod.get();
		modc->choked_ = true;
		add(mod);

		ModulePtr mod1(new chokedConv("Convergent", nodes_[0].get(), nodes_[1].get(), chokedConv::CALCULATED));
		add(mod1);

	}

	Inlet::Inlet() :
			Collection() {
	}

	GPObject &Inlet::duplicate() {
		Inlet *newMod(new Inlet(*this));
		duplicate((Collection *) newMod);

		return *newMod;
	}

	Inlet::Inlet(std::istream &is, const NodeMap &nodeMap, const ModuleMap &moduleMap)
			:
			Collection(is, nodeMap, moduleMap) {
	}

	Inlet::~Inlet() {
		// TODO Auto-generated destructor stub
	}

	void Inlet::reduce(const double &A) {
		N1_->A(A);
	}

	double Inlet::reduce() const {
		return N1_->A();
	}

	bool Inlet::isreduceable() const {
		return true;
	}

	void Inlet::unreduce() {
		N1_->A(N1_->Amax_);
	}

	double Inlet::drag() const {
		//here we should add also the non adaption drag
		//This definition is good for optimization but not for normal calculations
		return N1_->A() * N1_->rho() * std::pow(N1_->getU(), 2.0);
	}

	Glib::RefPtr<ModuleGraph> Inlet::draw() {
		graph_ = ModuleGraph::create(this);
		return graph_;
	}

	Inlet::OutType Inlet::out() const {
		OutType out(Collection::out());
		std::pair<std::string, std::vector<double> > tmp(out[1]);
		out[1] = out[2];
		out[2] = tmp;

		return out;
	}

bool Inlet::inputFreeStream()const{
	return true;
}

void Inlet::unserialize(const Archive& ar){
	Collection::unserialize(ar);

	//Make sure the input/ouput nodes are linked to this module and not
	//to its child modules
	N1(N1());
	N2(N2());
}
}
