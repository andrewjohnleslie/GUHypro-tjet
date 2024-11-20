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
#include "Rocket.h"
#include "solvers/isentropicDuct.h"

namespace hypro {
	Rocket::Rocket(std::string name, Node *N1, Node *N2, int ID)
			:
			propModule(name, N1, N2, ID, 0) {
	}

	Rocket::Rocket(std::istream &is, const NodeMap &nodeMap, const ModuleMap &moduleMap)
			:
			propModule(is, nodeMap, moduleMap) {
		try {
			is >> Isp_ >> Me_ >> mfr_ >> mfrMax_;
		} catch (std::exception &e) {
			throw std::runtime_error("Error: Rocket input incorrect\nProvide Isp_ Me_ mfr_ mfrMax_.");
		}
	}

	Rocket::Rocket() :
			propModule(),
			mfr_(0),
			mfrMax_(0),
			Isp_(0),
			Me_(0) {
	}

	Rocket::Rocket(const Rocket &mod)
			:
			propModule(mod),
			mfr_(mod.mfr_),
			mfrMax_(mod.mfrMax_),
			Isp_(mod.Isp_),
			Me_(mod.Me_) {
	}

	Rocket::~Rocket() {
		// TODO Auto-generated destructor stub
	}

	double Rocket::calculate() {
		//N2_->T_ = 300; // Only to asses gamma
		N2_->setTemp(300);
		N1_->A(0); //Useful in case of optimisation

		double p_p0;
		double T_T0;
		isentropicDuct::isentropic(Me_, N2_->gamma(), p_p0, T_T0);
		double pstr_p0;
		double Tstr_T0;
		isentropicDuct::isentropic(1, N2_->gamma(), pstr_p0, Tstr_T0);
		double T_Tstr = T_T0 / Tstr_T0;
		double Cf = std::sqrt(T_Tstr) * (Me_ + 1 / (N2_->gamma() * Me_));
		double cstr = Isp_ * g0_ / Cf;
		double Tstr = std::pow(cstr, 2.0) /
					  (N2_->gamma() * (N2_->R() / N2_->W())); //Specific Gas Constant i.e. R_spec = N2_->R() / N2_->W()
		//N2_->T_ = T_Tstr*Tstr;
		//N2_->p_ = 1e5; //Dummy value just to avoid error check
		N2_->setTP(T_Tstr * Tstr, 1e5);
		// N2_->U_ = Me_*N2_->geta();
		N2_->setU(Me_ * N2_->geta());
		//double rho = mfr_/(N2_->A()*N2_->U_);
		double rho = mfr_ / (N2_->A() * N2_->getU());
		//N2_->p_ = rho*N2_->R()*N2_->T_;
		N2_->setTP(T_Tstr * Tstr, rho * (N2_->R() / N2_->W()) * N2_->getTemp());

		return 1.0; //Never choked
	}

	std::vector<double> Rocket::propellantMfr() const {
		std::vector<double> Mfr(N2_->NfrX());

		std::vector<double> w(N2_->WX());
		std::vector<double> Coeff(Mfr.size(), 0.0);
		Coeff[0] = 1.0;
		Coeff[3] = 0.5;
		for (unsigned i = 0; i < Mfr.size(); i++) {
			Mfr[i] = Mfr[5] * Coeff[i] * w[i];
		}
		return Mfr;
	}

	void Rocket::reduce(const double &mfr) {
		if (mfr > mfrMax_) {
			throw std::runtime_error("Error: mfr cannot be higher than mfrMax");
		}
		if (mfr < 0) {
			throw std::runtime_error("Error: mfr cannot be lower than 0");
		}
		mfr_ = mfr;
	}

	double Rocket::reduce() const {
		return mfr_;
	}

	bool Rocket::isreduceable() const {
		return true;
	}

	void Rocket::unreduce() {
		mfr_ = mfrMax_;
	}

void Rocket::serialize(Archive& ar) const{
	propModule::serialize(ar);

	ar.put("Specific impulse", Isp_);
	ar.put("Exit Mach", Me_);
	ar.put("Mass flow rate", mfr_);
	ar.put("Maximum mass flow rate", mfrMax_);
}

void Rocket::unserialize(const Archive& ar) {
	propModule::unserialize(ar);

	Isp_ = ar.getDouble("Specific impulse");
	Me_ = ar.getDouble("Exit Mach");
	mfr_ = ar.getDouble("Mass flow rate");
	mfrMax_ = ar.getDouble("Maximum mass flow rate");
}

bool Rocket::canChoke()const{
	return false;
}

std::vector<const Node*> Rocket::inNode()const{
	return std::vector<const Node*>();
}
}


//BOOST_CLASS_EXPORT_IMPLEMENT(Rocket);
