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
#include "EffPhiComb.h"
#include "solvers/balanceMach.h"
#include "solvers/Friction.h"
#include "interpolation.h"
//#include <boost/serialization/vector.hpp>

namespace hypro {
    template<class Delta>
    EffPhiComb<Delta>::EffPhiComb(Node *N1, Node *N2)
            :
            EffPhiComb<Delta>(this->defName_, N1, N2) {
    }

    template<class Delta>
    EffPhiComb<Delta>::EffPhiComb(std::string name, Node *N1, Node *N2)
            :
            Combustion<Delta>(name, N1, N2, this->defID_),
            iFuel_(0) {
    }

    template<class Delta>
    EffPhiComb<Delta>::EffPhiComb(std::istream &is,
                                  const propModule::NodeMap &nodeMap, const propModule::ModuleMap &moduleMap)
            :
            Combustion<Delta>(is, nodeMap, moduleMap) {
        is >> iFuel_;
    }

    template<class Delta>
    EffPhiComb<Delta>::EffPhiComb::EffPhiComb()
            :
            Combustion<Delta>(),
            iFuel_(0),
            eta_() {
    }

    template<class Delta>
    EffPhiComb<Delta>::EffPhiComb::EffPhiComb(const EffPhiComb &mod)
            :
            Combustion<Delta>(mod) {
    }

    template<class Delta>
    EffPhiComb<Delta>::~EffPhiComb() {
        // TODO Auto-generated destructor stub
    }

    template<class Delta>
    void EffPhiComb<Delta>::delta(
            double &dG,
            double &dI,
            double &dH) const {
        Combustion<Delta>::delta(dG, dI, dH);

        double num = 0.0;
        double den = 0.0;
        int iOxider;
        for (std::size_t i = 0; i < this->Rcoeff_.size(); i++) {
            if (i == iFuel_) {
                num = this->Rcoeff_[i];
            } else if (this->Rcoeff_[i] != 0) {
                den = this->Rcoeff_[i];
                iOxider = i;
            }
        }
        if (num == 0.0) {
            throw std::logic_error("Error: fuel not present in the stoichiometric coefficients.");
        }
        if (den == 0.0) {
            throw std::logic_error("Error: oxider not present in the stoichiometric coefficients.");
        }
        double stoichComp = num / den;
        double phi = (this->N1().X()[iFuel_] / this->N1().X()[iOxider]) / stoichComp;

        double eta = interpolation::linear(eta_[0], eta_[1], phi);
        Node N1ref(*this->N1_);
        Node N2ref(*this->N2_);
        N1ref.setTemp(300);
        N2ref.setTemp(300);

        //dH += (1-eta)*this->N1_->G()*(this->N1_->gas().H(300) - this->N2_->gas().H(300));
        dH += (1 - eta) * this->N1_->G() * (N1ref.h() - N2ref.h()); //TODO Check if h or H
    }

template<class Delta>
void EffPhiComb<Delta>::serialize(Archive& ar) const{
	Combustion<Delta>::serialize(ar);

	ar.put("Fuel specie ID", iFuel_);
}

template<class Delta>
void EffPhiComb<Delta>::unserialize(const Archive& ar) {
	Combustion<Delta>::unserialize(ar);

	iFuel_ = ar.get<int>("Fuel specie ID");
}

template<>
std::string EffPhiComb<balanceMach>::typeName()const{
	return "EffPhiComb";
}

template<>
std::string EffPhiComb<Friction>::typeName()const{
	return "EffPhiCombFriction";
}

//Template implementation
    template
    class EffPhiComb<balanceMach>;


    template
    class EffPhiComb<Friction>;
}
//BOOST_CLASS_EXPORT_IMPLEMENT(EffPhiComb<balanceMach>);
//BOOST_CLASS_EXPORT_IMPLEMENT(EffPhiComb<Friction>);
