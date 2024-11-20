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
#include "Mixer.h"

namespace hypro {
	Mixer::Mixer(std::string name, Node *N1, Node *N2, Node *N3, int ID)
			:
			PropModule3(name, N1, N2, N3, ID),
			teta_(0.0),
			eff_(1.0){
}

	Mixer::Mixer() :
			PropModule3(),
			teta_(0.0),
			eff_(1.0) {
	}

	Mixer::Mixer(std::istream &is, const NodeMap &nodeMap, const ModuleMap &moduleMap)
			:
			PropModule3(is, nodeMap, moduleMap) {
		is >> eff_;
	}

	Mixer::Mixer(const Mixer &mod)
			:
			PropModule3(mod),
			teta_(mod.teta_),
			eff_(mod.eff_) {
	}

	Mixer::~Mixer() {
		// TODO Auto-generated destructor stub
	}

	void Mixer::changeComposition() {
		double N2 = N1_->Nfr() + N3_->Nfr();
		N2_->X((N1_->Nfr() * N1_->X() + N3_->Nfr() * N3_->X()) / N2);
	}

	void Mixer::delta(
			double &dG,
			double &dI,
			double &dH) const {
		dG = -N3_->A() / N1_->A() * N3_->rho() * N3_->getU();

		/* The hypothesis is that inlet 3 is close to station 1, so the
        pressure surrounding the inlet 3 is p1*/
		dI = eff_ * (dG * N3_->getU() - N3_->A() / N1_->A() * (N3_->getPress() - N1_->getPress())) *
			 std::cos(teta_ * pi_ / 180);

		dH = dG * (0.5 * std::pow(N3_->getU(), 2.0) + N3_->h());
	}

void Mixer::serialize(Archive& ar) const{
	PropModule3::serialize(ar);

	ar.put("Efficiency", eff_);
	ar.put("Theta", teta_);

}

void Mixer::unserialize(const Archive& ar) {
	PropModule3::unserialize(ar);

	eff_ = ar.getDouble("Efficiency");
}
}
