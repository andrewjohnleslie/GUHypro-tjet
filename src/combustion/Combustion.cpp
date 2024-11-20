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
#include "Combustion.h"
#include "solvers/balanceMach.h"
#include "solvers/Friction.h"

//#include <boost/serialization/vector.hpp>

namespace hypro {
template<class Delta>
Combustion<Delta>::Combustion(std::string name, Node* N1, Node* N2, int ID)
:
			Delta(name,N1,N2,ID,2),
			Rcoeff_(N1->X().size(),0.0),
			Pcoeff_(N1->X().size(),0.0){
}

	template<class Delta>
	Combustion<Delta>::Combustion():
			Delta(this->defName_, NULL, NULL, this->defID_, 2),
			Rcoeff_(),
			Pcoeff_() {
	}

	template<class Delta>
	Combustion<Delta>::Combustion(std::istream &is,
								  const propModule::NodeMap &nodeMap, const propModule::ModuleMap &moduleMap)
			:
			Delta(is, nodeMap, moduleMap),
			Rcoeff_(this->N1_->X().size(), 0),
			Pcoeff_(this->N1_->X().size(), 0) {
		try {
			for (std::size_t i = 0; i < Rcoeff_.size(); i++) {
				is >> Rcoeff_[i];
			}
			for (std::size_t i = 0; i < Pcoeff_.size(); i++) {
				is >> Pcoeff_[i];
			}
		} catch (std::exception &e) {
			throw std::runtime_error(
					"Error: Combustion Input Error. Requires Number of Rcoeff/Pcoeff equal to no. chemical species");
		}
	}


	template<class Delta>
	Combustion<Delta>::Combustion(const Combustion &mod)
			:
			Delta(mod),
			Rcoeff_(mod.Rcoeff_),
			Pcoeff_(mod.Pcoeff_) {
	}

	template<class Delta>
	Combustion<Delta>::~Combustion() {
		// TODO Auto-generated destructor stub
	}

	template<class Delta>
	void Combustion<Delta>::changeComposition() {
		//std::cout << this->N1_->report() << endl;
		if (Xswitch_ == 0)
		{

			std::vector<double>::const_iterator it_max = std::max_element(Rcoeff_.begin(), Rcoeff_.end());
			int i_max = it_max - Rcoeff_.begin();

			double limit = 1.0;
			double i_lim = 0;
			for (std::size_t i = 0; i < Rcoeff_.size(); i++) {
				if (Rcoeff_[i] != 0) {
					double limit1 = *it_max * this->N1_->X()[i] / Rcoeff_[i] * this->N1_->X()[i_max];
					if (limit > limit1) {
						limit = limit1;
						i_lim = i;
					}
				}
			}
			double Xlim = this->N1_->X()[i_lim] / Rcoeff_[i_lim]; //Normalized with its stochiometric coeff.

			std::vector<double> X(this->N1_->X().size());
			for (std::size_t i = 0; i < X.size(); i++) {
				X[i] = this->N1_->X()[i] - Xlim * Rcoeff_[i] + Xlim * Pcoeff_[i];
			}
			this->N2_->X(X);
			//std::cout << this->N2_->report() << endl;
		}
		if (Xswitch_ == 1)
		{

			this->N2_->X(XSET_);

		}
	}


template<class Delta>
void Combustion<Delta>::serialize(Archive& ar) const{
	Delta::serialize(ar);

	ar.put("Reactants coefficients", Rcoeff_);
	ar.put("Products coefficients", Pcoeff_);
}

template<class Delta>
void Combustion<Delta>::unserialize(const Archive& ar){
	Delta::unserialize(ar);

	Rcoeff_ = ar.getVector("Reactants coefficients");
	Pcoeff_ = ar.getVector("Products coefficients");
}


template<>
std::string Combustion<balanceMach>::typeName()const{
	return "Combustion";
}

template<>
std::string Combustion<Friction>::typeName() const {
		return "CombFriction";
	}

//Template implementation
template class Combustion<balanceMach>;

template class Combustion<Friction>;
}

//BOOST_CLASS_EXPORT_IMPLEMENT(Combustion<balanceMach>);
//BOOST_CLASS_EXPORT_IMPLEMENT(Combustion<Friction>);
