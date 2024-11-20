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
#include "InjectionPlate.h"
//#include <boost/serialization/vector.hpp>

namespace hypro {
	InjectionPlate::InjectionPlate(std::string name, Node *N1, Node *N2, int ID)
			:
			propModule(name, N1, N2, ID, 0),
			mfr_(0),
			mfrMax_(0),
			Y_(),
			T_(0),
			M_(0) {
	}

	InjectionPlate::InjectionPlate()
			:
			propModule(),
			mfr_(0),
			mfrMax_(0),
			Y_(),
			T_(0),
			M_(0) {
	}

	InjectionPlate::InjectionPlate(const InjectionPlate &mod)
			:
			propModule(mod),
			mfr_(mod.mfr_),
			mfrMax_(mod.mfrMax_),
			Y_(mod.Y_),
			T_(mod.T_),
			M_(mod.M_) {
	}

InjectionPlate::InjectionPlate(std::istream& is, const NodeMap& nodeMap, const ModuleMap& moduleMap)
:
		propModule(is,nodeMap,moduleMap),
		Y_(this->N1_->X().size(),0.0){
	try{
		is >> mfr_ >> mfrMax_;
		for(std::size_t i=0; i<Y_.size(); i++){
			is >> Y_[i];
		}
		is >> T_ >> M_;
	}catch(std::exception& e){
		throw std::runtime_error("Error: InjectionPlate input incorrect\nExpect mfr_ MfrMax_ injected_species temperature mach number");
	}
}

	InjectionPlate::~InjectionPlate() {
		// TODO Auto-generated destructor stub
	}

	double InjectionPlate::calculate() {

		N1_->A(0); //Useful in case of optimisation

		N2_->setTPY(T_, 10, Y_);
		N2_->setU(M_ * N2_->geta());
		double rho = mfr_ / (N2_->getU() * N2_->A());
		try {
			N2_->setPress(rho * N2_->RT());
		} catch (Cantera::CanteraError e) {
			cout << "Density zero" << endl;
		}
		//TODO check if setTemp and setPress actually work

		return 1.0; //Never choked
	}

	std::vector<double> InjectionPlate::propellantMfr() const {
		return N2_->mfrX();
	}

	void InjectionPlate::reduce(const double &mfr) {
		if (mfr > mfrMax_) {
			throw std::runtime_error("Error: mfr cannot be higher than mfrMax");
		}
		if (mfr < 0) {
			throw std::runtime_error("Error: mfr cannot be lower than 0");
		}
		mfr_ = mfr;
	}

	double InjectionPlate::reduce() const {
		return mfr_;
	}

	bool InjectionPlate::isreduceable() const {
		return true;
	}

	void InjectionPlate::unreduce() {
		mfr_ = mfrMax_;
	}

void InjectionPlate::serialize(Archive& ar) const{
	propModule::serialize(ar);

	ar.put("Mass flow rate", mfr_);
	ar.put("Maximum mass flow rate", mfrMax_);
	ar.put("Composition of injected mixture", Y_);
	ar.put("Temperature of injection", T_);
	ar.put("Mach", M_);
}

void InjectionPlate::unserialize(const Archive& ar) {
	propModule::unserialize(ar);

	mfr_ = ar.getDouble("Mass flow rate");
	mfrMax_ = ar.getDouble("Maximum mass flow rate");
	Y_ = ar.getVector("Composition of injected mixture");
	T_ = ar.getDouble("Temperature of injection");
	M_ = ar.getDouble("Mach");
}

bool InjectionPlate::canChoke()const{
	return false;
}

std::vector<const Node*> InjectionPlate::inNode()const{
	return std::vector<const Node*>();
}
}
