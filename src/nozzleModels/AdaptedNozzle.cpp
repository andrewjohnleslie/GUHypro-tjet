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
#include "AdaptedNozzle.h"
#include <solvers/EffIsenDuct.h>

namespace hypro {
	template<class Eff>
	AdaptedNozzle<Eff>::AdaptedNozzle(std::string name, Node *N1, Node *N2)
			:
			AdaptedNozzle(name, N1, N2, this->defID_) {
	}

	template<class Eff>
	AdaptedNozzle<Eff>::AdaptedNozzle(std::string name, Node *N1, Node *N2, int ID)
			:
			Eff(name, N1, N2, ID),
			freeStream_(0) {
	}

	template<class Eff>
	AdaptedNozzle<Eff>::AdaptedNozzle():
			Eff(),
			freeStream_(0) {
	}

	template<class Eff>
	AdaptedNozzle<Eff>::AdaptedNozzle(std::istream &is,
									  const propModule::NodeMap &nodeMap, const propModule::ModuleMap &moduleMap)
			:
			Eff(is, nodeMap, moduleMap) {
		try {
			unsigned FSid;
			is >> FSid;
			freeStream_ = nodeMap.at(FSid).get();
		} catch (std::exception &e) {
			throw std::runtime_error("Error: Adapted Nozzle Input Incorrect");
		}
	}

	template<class Eff>
	AdaptedNozzle<Eff>::AdaptedNozzle(const AdaptedNozzle &mod)
			:
			Eff(mod),
			freeStream_(mod.freeStream_) {
	}

	template<class Eff>
	AdaptedNozzle<Eff>::~AdaptedNozzle() {
		// TODO Auto-generated destructor stub
	}

	template<class Eff>
	double AdaptedNozzle<Eff>::calculate() {
		double p1_p01, T1_T01;

		this->choked_ = !(this->N1_->M() > 1.0);

		this->isentropic(this->N1_->M(), this->N1_->gamma(), p1_p01, T1_T01);

		this->pressLoss();
		this->thermLoss();
		this->applyLoss();

		this->N2_->setPress(freeStream_->getPress());
		this->N2_->X(this->N1_->X());

		/**
            Outside try loop is to take care of really low pressures, which cause issues for isentropicP.
            The catch tests if the pressure is low, and if so does not try and adapt the nozzle.
        **/
		this->N2_->isentropicP(this->N1mod(), this->N2_->getPress());
		double dH = this->N1mod().H() - this->N2_->h();
		if (dH >= 0.0) {
			this->N2_->setU(std::sqrt(2 * dH));
			try {
				double M2 = this->N2_->M();
				this->N2_->A(this->N1_->A() * this->mfrEq(M2));
				return 1; //not choked
			} catch (const std::range_error &e) {
				this->N2_->A(this->N2_->Amax_);
				return isentropicDuct::calculate();
			}
		} else { //In case the atmospheric values are higher than the stagnation just open at maximum
			this->N2_->A(this->N2_->Amax_);
			return isentropicDuct::calculate();
		}
	}


template<>
std::string AdaptedNozzle<isentropicDuct>::typeName() const{
	return "AdaptedNozzle";
}

	template<>
	std::string AdaptedNozzle<EffIsenDuct>::typeName() const {
		return "AdaptedNozzEff";
	}


template<class Eff>
void AdaptedNozzle<Eff>::serialize(Archive& ar) const{
	Eff::serialize(ar);

	ar.put("Free Stream node", freeStream_->ID_);
}

template<class Eff>
void AdaptedNozzle<Eff>::unserialize(const Archive& ar) {
	Eff::unserialize(ar);

	freeStream_ = ar.getNodeRef("Free Stream node");
}

//Template implementation
template class AdaptedNozzle< isentropicDuct >;
template class AdaptedNozzle< EffIsenDuct >;
}
