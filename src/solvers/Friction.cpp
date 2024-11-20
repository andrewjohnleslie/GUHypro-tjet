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
#include "Friction.h"

namespace hypro {
	Friction::Friction(std::string name, Node *N1, Node *N2, int ID)
			:
			balanceMach(name, N1, N2, ID, 2),
			Cf_(0.0),
			L_(0.0),
			D_(0.0){
}

Friction::Friction(std::string name, Node *N1, Node *N2, int ID, int args)
			:
			balanceMach(name, N1, N2, ID, args),
			Cf_(0.0),
			L_(0.0),
			D_(0.0){
			}

Friction::Friction():
	balanceMach(),
	Cf_(0.0),
	L_(0.0),
	D_(0.0){
}

Friction::Friction(std::istream& is, const NodeMap& nodeMap, const ModuleMap& moduleMap)
:
		balanceMach(is,nodeMap,moduleMap){
	is >> Cf_ >> L_ >> D_;
}


	Friction::Friction(const Friction &mod)
			:
			balanceMach(mod),
			Cf_(mod.Cf_),
			L_(mod.L_),
			D_(mod.D_) {
	}

	Friction::~Friction() {
		// TODO Auto-generated destructor stub
	}

	void Friction::changeComposition() {
		N2_->X(N1_->X());
	}

	void Friction::delta(
			double &dG,
			double &dI,
			double &dH) const {
		dG = 0.0;
		dI = 2 * Cf_ * (L_ / D_) *
			 ((N1_->rho() * std::pow(N1_->getU(), 2.0)) /*+ (N2_->rho()*std::pow(N2_->U_,2.0)) */);
		if (heat_flux){

			dH = (qdot)/(N1_->A()*2);
			// std::cout << "dH: "   << dH << std::endl;
		} else {
			dH = 0.0;
		}

	}

	std::string Friction::typeName() const {

		std::string name = "Friction";
		return name;
	}

	void Friction::serialize(Archive& ar) const{
	balanceMach::serialize(ar);

	ar.put("Friction coefficient", Cf_);
	ar.put("Length", L_);
	ar.put("Hydraulic diameter", D_);
}

void Friction::unserialize(const Archive& ar) {
	balanceMach::unserialize(ar);

	Cf_ = ar.getDouble("Friction coefficient");
	L_ = ar.getDouble("Length");
	D_ = ar.getDouble("Hydraulic diameter");
}
}
