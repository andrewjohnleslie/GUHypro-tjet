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

#include "utilities/Archive.h"
#include "ParMixer.h"

namespace hypro {
	ParMixer::ParMixer(std::string name, Node *N1, Node *N2, Node *N3, int ID)
			:
			PropModule3(name, N1, N2, N3, ID),
			teta_(0.0),
			eff_(1.0){
}

	ParMixer::ParMixer()
			:
			PropModule3(),
			teta_(0.0),
			eff_(1.0) {
	}

	ParMixer::ParMixer(std::istream &is, const NodeMap &nodeMap, const ModuleMap &moduleMap)
			:
			PropModule3(is, nodeMap, moduleMap) {
		is >> eff_;
	}

	ParMixer::ParMixer(const ParMixer &mod)
			:
			PropModule3(mod),
			teta_(mod.teta_),
			eff_(mod.eff_) {
	}

	ParMixer::~ParMixer() {
		// TODO Auto-generated destructor stub
	}

	void ParMixer::checkAreas() {
		N2_->A(N1_->A() + N3_->A());
		N2_->Amax_ = N2_->A();
		N2_->Amin_ = N2_->A();
	}

	void ParMixer::changeComposition() {
		double N2 = N1_->Nfr() + N3_->Nfr();
		N2_->X((N1_->Nfr() * N1_->X() + N3_->Nfr() * N3_->X()) / N2);
	}

	void ParMixer::delta(
			double &dG,
			double &dI,
			double &dH) const {
//		std::cout << N1_->mfr() << std::endl;
//		std::cout << N3_->mfr() << std::endl;
		dG = -N3_->A() / N2_->A() * N3_->rho() * N3_->getU();

		dI = eff_ * (dG * N3_->getU() - N3_->A() / N2_->A() * N3_->getPress()) * std::cos(teta_ * pi_ / 180);


	//Efficiency applied to node 1 too
	double dG1 = -N1_->A()/N2_->A()*N1_->rho()*N1_->getU();
	dI -= (1 - eff_)*(dG1*N1_->getU() - N1_->A()/N2_->A()*N1_->getPress());

	dH = dG*(0.5*std::pow(N3_->getU(),2.0) + N3_->h());
}

void ParMixer::serialize(Archive& ar) const{
	PropModule3::serialize(ar);

	ar.put("Efficiency", eff_);
}

void ParMixer::unserialize(const Archive& ar) {
	PropModule3::unserialize(ar);

	eff_ = ar.getDouble("Efficiency");
}
}
