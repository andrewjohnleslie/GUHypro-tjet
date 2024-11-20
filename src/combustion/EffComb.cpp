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
#include "EffComb.h"
#include "solvers/balanceMach.h"
#include "solvers/Friction.h"


namespace hypro {
	template<class Delta>
	EffComb<Delta>::EffComb(Node *N1, Node *N2)
			:
			EffComb<Delta>(this->defName_, N1, N2) {
	}

	template<class Delta>
	EffComb<Delta>::EffComb(std::string name, Node *N1, Node *N2)
			:
			Combustion<Delta>(name, N1, N2, this->defID_),
			eta_(1.0) {
	}

	template<class Delta>
	EffComb<Delta>::EffComb():
			Combustion<Delta>(),
			eta_(1.0) {
	}

	template<class Delta>
	EffComb<Delta>::EffComb(std::istream &is,
							const propModule::NodeMap &nodeMap, const propModule::ModuleMap &moduleMap)
			:
			Combustion<Delta>(is, nodeMap, moduleMap) {
		is >> eta_;
	}

	template<class Delta>
	EffComb<Delta>::EffComb::EffComb(const EffComb &mod)
			:
			Combustion<Delta>(mod),
			eta_(mod.eta_) {
	}

	template<class Delta>
	EffComb<Delta>::~EffComb() {
		// TODO Auto-generated destructor stub
	}

	template<class Delta>
	void EffComb<Delta>::delta(
			double &dG,
			double &dI,
			double &dH) const {
		Combustion<Delta>::delta(dG, dI, dH);

		Node N1ref(*this->N1_);
		Node N2ref(*this->N2_);
		N1ref.setTemp(300);
		N2ref.setTemp(300);

		//dH += (1-eta_)*this->N1_->G()*(this->N1_->gas().H(300) - this->N2_->gas().H(300));
		dH += (1 - eta_) * this->N1_->G() * (N1ref.h() - N2ref.h()); //TODO Check if h or H
	}

template<class Delta>
void EffComb<Delta>::serialize(Archive& ar) const{
	Combustion<Delta>::serialize(ar);

	ar.put("Efficiency", eta_);
}

template<class Delta>
void EffComb<Delta>::unserialize(const Archive& ar) {
	Combustion<Delta>::unserialize(ar);

	eta_ = ar.getDouble("Efficiency");
}

template<>
std::string EffComb<balanceMach>::typeName()const{
	return "EffComb";
}

template<>
std::string EffComb<Friction>::typeName()const{
	return "EffCombFriction";
}

//Template implementation
	template
	class EffComb<balanceMach>;


	template
	class EffComb<Friction>;
}
//BOOST_CLASS_EXPORT_IMPLEMENT(EffComb<balanceMach>);
//BOOST_CLASS_EXPORT_IMPLEMENT(EffComb<Friction>);

