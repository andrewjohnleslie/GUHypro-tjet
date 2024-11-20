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
#include "InjectionPhiNoTemp.h"
#include "solvers/balanceMach.h"
#include "solvers/Friction.h"
#include <InfNanFloat.h>
//#include <boost/serialization/vector.hpp>
namespace hypro {

	template<class Delta>
	InjectionPhiNoTemp<Delta>::InjectionPhiNoTemp(std::string name, Node *N1, Node *N2)
			:
			InjectionPhiNoTemp<Delta>(name, N1, N2, this->defID_) {
	}

template<class Delta>
InjectionPhiNoTemp<Delta>::InjectionPhiNoTemp(std::string name, Node* N1, Node* N2, int ID)
:
		Delta(name,N1,N2,ID,2),
		Rcoeff_(N1->X().size(),0),
		Phi_(1.0),
		PhiMax_(1.0),
		mfrMax_(std::numeric_limits<double>::infinity()),
		Xf_(N1->X().size(),0.0),
		Tf_(300.0){
}

	template<class Delta>
	InjectionPhiNoTemp<Delta>::InjectionPhiNoTemp()
			:
			Delta(),
			Rcoeff_(),
			Phi_(1.0),
			PhiMax_(1.0),
			mfrMax_(std::numeric_limits<double>::infinity()),
			Xf_(),
			Tf_(300.0) {
	}

	template<class Delta>
	InjectionPhiNoTemp<Delta>::InjectionPhiNoTemp(std::istream &is,
									  const propModule::NodeMap &nodeMap, const propModule::ModuleMap &moduleMap)
			:
			Delta(is, nodeMap, moduleMap),
			Rcoeff_(this->N1_->X().size(), 0),
			Xf_(this->N1_->X().size(), 0.0) {
		try {
			for (std::size_t i = 0; i < Rcoeff_.size(); i++) {
				is >> Rcoeff_[i];
			}
			is >> Phi_ >> PhiMax_ >> mfrMax_;
			for (std::size_t i = 0; i < Xf_.size(); i++) {
				is >> Xf_[i];
			}
			is >> Tf_;
		} catch (std::exception &e) {
			throw std::runtime_error("Error: Incorrect inputs.\nProvide Rcoeff, Phi, PhiMax, mfrMax, Xf_, Tf_");
		}
//}catch(std::exception& e){
//		throw std::runtime_error("Error: Rocket input incorrect\nProvide Isp_ Me_ mfr_ mfrMax_.");
//	}
	}

	template<class Delta>
	InjectionPhiNoTemp<Delta>::InjectionPhiNoTemp(const InjectionPhiNoTemp &mod)
			:
			Delta(mod),
			Rcoeff_(mod.Rcoeff_),
			Phi_(mod.Phi_),
			PhiMax_(mod.PhiMax_),
			mfrMax_(mod.mfrMax_),
			Xf_(mod.Xf_),
			Tf_(mod.Tf_) {
	}

	template<class Delta>
	InjectionPhiNoTemp<Delta>::~InjectionPhiNoTemp() {
		// TODO Auto-generated destructor stub
	}

	template<class Delta>
	double InjectionPhiNoTemp<Delta>::findEqv() const {

		double denom;
		double num;

		for (std::size_t i = 0; i < Xf_.size(); i++) {
			if (Xf_[i] != 0) {
				num = Rcoeff_[i];

			}
			if (Xf_[i] == 0 && Rcoeff_[i] != 0) {
				denom = Rcoeff_[i];
			}
		}
		double eqv = num / denom;
		return eqv;
	}

	template<class Delta>
	void InjectionPhiNoTemp<Delta>::changeComposition() {
		double stoichcomp = findEqv();
		double K = (stoichcomp) * Phi_ * this->N1_->X("O2");  //2 = h-o
		vector<double> Xfupdate = Xf_ * K;
		//vector<double> Xinitial = this->N1_->X();
		this->N2_->X((this->N1_->X() + Xfupdate) / (1 + K)); //TODO Transform for multiplying contents of two vectors
	}

	template<class Delta>
	void InjectionPhiNoTemp<Delta>::delta(
			double &dG,
			double &dI,
			double &dH) const {
		Delta::delta(dG, dI, dH);

		Node Nfuel(*this->N1_);
		Nfuel.setTPX(Tf_, Nfuel.getPress(), Xf_);
		//double Wf = this->N1_->specieGas().getMixture(Xf_).W();
		double Wf = Nfuel.W();
		double stoichcomp = findEqv();
		dG += -(stoichcomp) * Phi_ * (Wf / this->N1_->W()) * this->N1_->X("O2") * this->N1_->rho() * this->N1_->getU();
//		std::cout << "InjectionPhi<Delta>::delta(): " << Phi_ << std::endl;
		dI += 0;
//		std::cout << "stoichcomp is: " << dG << "\n";


		double hf = Nfuel.h();
		//specieThermo< janafThermo<perfectGas> > gas = this->N1_->specieGas().getMixture(Xf_);

		/*if(Tf_<gas.Tlow()){
            hf = gas.H(gas.Tlow()) + gas.Cp(gas.Tlow())*(Tf_ - gas.Tlow());
        }else if(Tf_<=gas.Thigh()){
            hf = gas.H(Tf_);
        }else{
            hf = gas.H(gas.Thigh()) + gas.Cp(gas.Thigh())*(Tf_ - gas.Thigh());
        }*/
//		std::cout << "dG is: " << dG << "\nhf is: " << hf << "\n";
		dH += 0 ;
	}

	template<class Delta>
	void InjectionPhiNoTemp<Delta>::reduce(const double &Phi) {
		if (Phi > PhiMax_) {
			//Foam::error("Error: Phi cannot be higher than PhiMax").abort();
			cout << "Error: Phi cannot be higher than PhiMax" << endl;
		}
		Phi_ = Phi;
//        std::cout << "InjectionPhi<Delta>::calculate(): " << Phi_ << std::endl;
	}

	template<class Delta>
	double InjectionPhiNoTemp<Delta>::reduce() const {
		return Phi_;
	}

	template<class Delta>
	bool InjectionPhiNoTemp<Delta>::isreduceable() const {
		return true;
	}


	template<class Delta>
	void InjectionPhiNoTemp<Delta>::unreduce() {
		Phi_ = PhiMax_;
//        std::cout << "InjectionPhi<Delta>::calculate(): " << Phi_ << std::endl;
	}

	template<class Delta>
	std::vector<double> InjectionPhiNoTemp<Delta>::propellantMfr() const {
		double dG, dummy;
		delta(dG, dummy, dummy);
		double mfr = -dG * this->N1_->A();
//    std::cout << "mfr = " << mfr << ", -dG: " << -dG << ", this_>N1_->A(): " << this->N1_->A() << std::endl;
		std::vector<double> w(this->N1_->WX());
		double W = 0.0;
		for (std::size_t i = 0; i < Xf_.size(); i++) {
			W = W + Xf_[i] * w[i];
		}
		double Nfr = mfr / W;

		std::vector<double> mfrX(Xf_);
		for (std::size_t i = 0; i < Xf_.size(); i++) {
			mfrX[i] = Xf_[i] * w[i] * Nfr;
		}

//    std::cout << "InjectionPhi<Delta>::propellantMfr(): ";
//    for (auto &massFlow : mfrX) {
//            std::cout << massFlow << ", ";
//    }
//    std::cout << "\n";
		return mfrX;
	}


template<class Delta>
void InjectionPhiNoTemp<Delta>::serialize(Archive& ar) const{
	Delta::serialize(ar);

	ar.put("Equivalence ratio (Phi)", Phi_);
	ar.put("Maximum Phi", PhiMax_);
	ar.put("Fuel composition", Xf_);
	ar.put("Reactants coefficients", Rcoeff_);
	ar.put("Temperature of injected fuel", Tf_);
	ar.put("Mass Flow Rate of fuel components", this->propellantMfr());
}

template<class Delta>
void InjectionPhiNoTemp<Delta>::unserialize(const Archive& ar) {
	Delta::unserialize(ar);

	Phi_ = ar.getDouble("Equivalence ratio (Phi)");
	PhiMax_ = ar.getDouble("Maximum Phi");
	Xf_ = ar.getVector("Fuel composition");
	Rcoeff_ = ar.getVector("Reactants coefficients");
	Tf_ = ar.getDouble("Temperature of injected fuel");
}

	template<>
	std::string InjectionPhiNoTemp<balanceMach>::typeName() const {
		return "InjectionPhiNoTemp";
	}

	template<>
	std::string InjectionPhiNoTemp<Friction>::typeName() const {
		return "InjPhiNoTempFriction";
	}

	template<class Delta>
	double InjectionPhiNoTemp<Delta>::calculate() {
		double dG, dummy;
		delta(dG, dummy, dummy);
		double mfr = -dG * this->N1_->A();
		Node Nfuel(*this->N1_);
		Nfuel.setTPX(Tf_, Nfuel.getPress(), Xf_); //TODO find out what this pressure should be
		if (mfr > mfrMax_) {
			//double Wf = this->N1_->specieGas().getMixture(Xf_).W();
			double Wf = Nfuel.W();
			double stoichcomp = findEqv();
			Phi_ = mfrMax_ / (stoichcomp * (Wf / this->N1_->W()) * this->N1_->X("O2") * this->N1_->mfr());
//            std::cout << "InjectionPhi<Delta>::calculate() Phi: " << Phi_ << std::endl;
		}

		return Delta::calculate();
	}


//Template implementation
	template
	class InjectionPhiNoTemp<balanceMach>;

	template
	class InjectionPhiNoTemp<Friction>;
}
