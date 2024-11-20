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
#include "InjectionPlatePressure.h"

namespace hypro {
    InjectionPlatePressure::InjectionPlatePressure(std::string name, Node *N1, Node *N2, int ID)
            :
            propModule(name, N1, N2, ID, 0),
            mfr_(0),
            mfrMax_(0),
            Y_(),
            T_(0),
            p_(0) {
    }

    InjectionPlatePressure::InjectionPlatePressure(const InjectionPlatePressure &mod)
            :
            propModule(mod),
            mfr_(mod.mfr_),
            mfrMax_(mod.mfrMax_),
            Y_(mod.Y_),
            T_(mod.T_),
            p_(mod.p_) {
    }

    InjectionPlatePressure::InjectionPlatePressure(std::istream &is, const NodeMap &nodeMap, const ModuleMap &moduleMap)
            :
            propModule(is, nodeMap, moduleMap),
            Y_(this->N1_->X().size(), 0.0) {
        try {
            is >> mfr_ >> mfrMax_;
            for (unsigned i = 0; i < Y_.size(); i++) {
                is >> Y_[i];
            }
            is >> T_ >> p_;
        } catch (std::exception &e) {
            throw std::runtime_error(
                    "Error: InjectionPlatePressure input incorrect\nExpect mfr_ MfrMax_ injected_species temperature mach number");
        }
    }

    InjectionPlatePressure::~InjectionPlatePressure() {
        // TODO Auto-generated destructor stub
    }

    double InjectionPlatePressure::calculate() {

        N1_->A(0); //Useful in case of optimisation
        N2_->setTPY(T_, p_, Y_);
        //N2_->setU(mfr_/((p_/N2_->RT())*N2_->A()));
        N2_->setU(mfr_ / (N2_->rho() * N2_->A()));
        // std::cout << mfr_ << std::endl;
//    std::cout << "rho man: " << p_/N2_->RT() << std::endl;
//    std::cout << "rho: " << N2_->rho() << std::endl;
//    std::cout << "U: " << mfr_/((p_/N2_->RT())*N2_->A()) << endl;
//    std::cout << "Manual: " << N2_->rho()*N2_->getU()*N2_->A() << std::endl;
//	std::cout << "mfr_: " << mfr_ << std::endl;
//	//N2_->setU(mfr_/(N2_->rho()*N2_->A()));
//    std::cout << "N2_->mfr():" << N2_->mfr() << std::endl;

        return 1.0; //Never choked
    }

    vector<double> InjectionPlatePressure::propellantMfr() const {
        return N2_->mfrX();
    }

    void InjectionPlatePressure::reduce(const double &mfr) {
        if (mfr > mfrMax_) {
            throw std::runtime_error("Error: mfr cannot be higher than mfrMax");
        }
        if (mfr < 0) {
            throw std::runtime_error("Error: mfr cannot be lower than 0");
        }
        mfr_ = mfr;
    }

    double InjectionPlatePressure::reduce() const {
        return mfr_;
    }

    bool InjectionPlatePressure::isreduceable() const {
        return true;
    }

    void InjectionPlatePressure::unreduce() {
        mfr_ = mfrMax_;
    }

    void InjectionPlatePressure::serialize(Archive& ar) const{
	propModule::serialize(ar);

	ar.put("Mass flow rate", mfr_);
	ar.put("Maximum mass flow rate", mfrMax_);
	ar.put("Composition of injected mixture", Y_);
	ar.put("Temperature of injection", T_);
	ar.put("Pressure", p_);
}

void InjectionPlatePressure::unserialize(const Archive& ar) {
	propModule::unserialize(ar);

	mfr_ = ar.getDouble("Mass flow rate");
	mfrMax_ = ar.getDouble("Maximum mass flow rate");
	Y_ = ar.getVector("Composition of injected mixture");
	T_ = ar.getDouble("Temperature of injection");
	p_ = ar.getDouble("Pressure");
}

bool InjectionPlatePressure::canChoke()const{
	return false;
}

std::vector<const Node*> InjectionPlatePressure::inNode()const{
	return std::vector<const Node*>();
}
}
