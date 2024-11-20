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
#include "isentropicDuct.h"
#include "nonLinear.h"

namespace hypro {
	isentropicDuct::isentropicDuct(std::string name, Node *N1, Node *N2, int ID)
			:
			propModule(name, N1, N2, ID, 2),
			N1mod_(),
			p02_p01_(1.0),
			T02_T01_(1.0),
			choked_(false){
}

isentropicDuct::isentropicDuct(std::istream& is, const NodeMap& nodeMap, const ModuleMap& moduleMap)
:
		propModule(is,nodeMap,moduleMap){
	try{
		is >> choked_;
	}catch(std::exception& e){
		throw std::runtime_error("Error: Incorrect Inputs");
	}
		}


	isentropicDuct::isentropicDuct(const isentropicDuct &mod)
			:
			propModule(mod),
			N1mod_(),
			p02_p01_(mod.p02_p01_),
			T02_T01_(mod.T02_T01_),
			choked_(mod.choked_) {
	}

	isentropicDuct::~isentropicDuct() {
		// TODO Auto-generated destructor stub
	}

	void isentropicDuct::pressLoss() {
		p02_p01_ = 1;
	}

	void isentropicDuct::thermLoss() {
		T02_T01_ = 1;

	}

	void isentropicDuct::applyLoss() {
		//TODO check with Alessandro how this is supposed to work - currently T02_T01_ is always 1 so it should work as it is now
		//Apply pressure drop
		thermoKineticState N1mod(*N1_);
		N1mod.setPress(N1_->getPress() * p02_p01_);

		//Apply temperature drop
		N1mod_.reset(new thermoKineticState(N1mod));
		thermoKineticState N1tmp(N1mod);
		N1tmp.setTemp(10); //Used to get the enthalpy at 0 Kelvin
		N1tmp.setTemp(298.15);
		// Modify T02_T01_ for our reference temperature: (T0_1 * ( T0_2/T0_1 ) - Tref)/(T0_1 = Tref)
		double T02_T01_mod = ((N1_->T0() * T02_T01_) - 298.15)/(N1_->T0() - 298.15);
		N1mod_->Th((N1mod.H()-N1tmp.h())*T02_T01_mod + N1tmp.h() - 0.5*std::pow(N1mod_->getU(),2.0));
		N1mod_->isentropicT(N1mod, N1mod_->getTemp());
	}

	const thermoKineticState &isentropicDuct::N1mod() const {
		return *N1mod_;
	}

	double isentropicDuct::calculate() {
		N2_->X(N1_->X());

		pressLoss();
		thermLoss();
		applyLoss();


        // Find throat area corresponding to chokingMach
		double M2 = chokingMach_;		// Value appears to be linked with turbojet.h value
		//std::cout << "IsenDuct: chokingMach (M2) = " << M2 << std::endl;	// 0.9 [Minf = 0.91] - when turbojet.h limit is set tpo that value
		double Astr = N1_->A() * mfrEq(M2);
		//std::cout << "IsenDuct: Astr = " << Astr << std::endl;	// 5.73486
		//std::cout << "IsenDuct: mfrEq() = " << mfrEq(M2) << std::endl;	// 1.00172
		double rChok = N2_->A() - Astr;
		//std::cout << "IsenDuct: rChoke = " << rChok << std::endl;	// 
		// if chokingMach throat area is greater than the actual throat, reduce the mass flow
		if (rChok < 0.0) {
//		    fmt::print("IsentropicDuct with Node 2: {} over-choked, reducing mass flow. r = {}\n", N2_->Name_, rChok);
			return rChok;
		}

		double Mmax, Mmin;
		if (choked_) {
			if (N1_->M() < chokingMach_) {
				Mmax = 10;
				Mmin = chokingMach_;
			} else {
				Mmax = chokingMach_;
				Mmin = 0;
			}
		} else {
			if (N1_->M() < chokingMach_) {
				Mmax = chokingMach_;
				Mmin = 0;
			} else {
				Mmax = 10;
				Mmin = chokingMach_;
			}
		}


//		M2 = nonLinear<isentropicDuct>::funSolve(&isentropicDuct::mfrEq, Mmin, Mmax, *this, 0);
//		mfrEq(M2);

        const boost::uintmax_t maxit = 30;
        boost::uintmax_t it = maxit;
        bool is_rising = true;
        int digits = std::numeric_limits<double>::digits;
        int get_digits = digits - 3;
        boost::math::tools::eps_tolerance<double> tol(get_digits);
        auto g_err = [this](double M2) {
            return mfrEq(M2, NULL);
        };

        std::pair<double, double> r = boost::math::tools::toms748_solve(g_err, Mmin, Mmax, tol,
                                                                                 it);
        M2 = r.first + (r.second - r.first) / 2;
        mfrEq(M2);

		if(verbosity_){
			fmt::print(N1_->report());
			fmt::print(N2_->report());
		}
		return rChok;
	}

	void isentropicDuct::isentropic(
			const double &Mach, const double &gamma,
			double &p_p0, double &T_T0) {
		T_T0 = 1 / (1 + 0.5 * (gamma - 1) * std::pow(Mach, 2.0));
		p_p0 = std::pow(T_T0, gamma / (gamma - 1));

	}

	double isentropicDuct::isentropicInv(
			const double &p_p0, const double &gamma) {
		if (p_p0 > 1.0) {
			throw std::runtime_error("Error: static pressure higher than total pressure.");
		}
		double T_T0 = std::pow(p_p0, (gamma - 1) / gamma);
		return std::sqrt((2 / (gamma - 1)) * (1 / T_T0 - 1.0));
	}

	double isentropicDuct::mfrEq(double &M2, const void *par[]) {
		if(verbosity_) {
			fmt::print("M2: {}\n", M2);
		}
		double mfreq = mfrEq(M2);
		if(verbosity_) {
			fmt::print("N2_->A(): {}, N1_->A(): {}, mfreq: {}, SUM: {} \n",  N2_->A() , N1_->A(), mfreq, N2_->A() / N1_->A() - mfreq);
		}
		return N2_->A() / N1_->A() - mfreq;
	}

	double isentropicDuct::mfrEqP(double &p2, const void *par[]) {

		if (p2 < toll_) { //To avoid fzero detect a solution near zero
			p2 = toll_;
		}

		N2_->isentropicP(*N1_, p2);
		N2_->setU(std::sqrt(2.0 * N1_->H() - N2_->h()));
                     
		return N2_->mfr() - N1_->mfr();
	}

	double isentropicDuct::pThroat(double &p, const void *par[]) {
		Node Nt(*N2_);
		Nt.isentropicP(*N1_, p);
		Nt.setU(std::sqrt(2.0 * N1_->H() - Nt.h()));
		return Nt.getU() - Nt.geta();
	}

	double isentropicDuct::mfrEq(double &M2) const {

		if (M2 < toll_) { //To avoid fzero detect geta solution near zero
			M2 = toll_;
		}

		if (!N1mod_.get()) {
			throw std::logic_error("Error: N1 modified state not assigned.");
		}

		N2_->setStagnation(N1mod_->getstag(), M2);
//
		if(verbosity_) {
			fmt::print("N1_->h(): {}, N2_->h(): {}, N1mod_->h(): {}\n", N1_->h(), N2_->h(), N1mod_->h());
			fmt::print("N1_->H(): {}, N2_->H(): {}, N1mod_->H(): {}\n", N1_->H(), N2_->H(), N1mod_->H());
			fmt::print("N1_->s(): {}, N2_->s(): {}, N1mod_->s(): {}\n", N1_->s(), N2_->s(), N1mod_->s());
			fmt::print("N1_->getU(): {}, N1_->rho: {}\n", N1_->getU(), N1_->rho());
			fmt::print("G1: {}, G2: {}, diffG: {}\n", N1_->G(), N2_->G(), N1_->G() / N2_->G());
		}


		return N1_->G() / N2_->G();
	}

void isentropicDuct::serialize(Archive& ar) const{
	propModule::serialize(ar);

	ar.put("Is choked", choked_);
}

void isentropicDuct::unserialize(const Archive& ar) {
	propModule::unserialize(ar);

	choked_ = ar.get<bool>("Is choked");
}
}
