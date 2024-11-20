/*!
 *  \author     Mark De Luca
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
#include "Fanno.h"
#include "nonLinear.h"

namespace hypro {
	Fanno::Fanno(Node *N1, Node *N2)
			:
			Fanno(defName_, N1, N2) {
	}

	Fanno::Fanno(std::string name, Node *N1, Node *N2)
			:
			propModule(name, N1, N2),
			Cf_(0.0),
			L_(0.0),
			D_(0.0) {
	}

	Fanno::Fanno(std::istream &is, const NodeMap &nodeMap, const ModuleMap &moduleMap)
			:
			propModule(is, nodeMap, moduleMap) {
		is >> Cf_ >> L_ >> D_;
	}

	Fanno::Fanno(const Fanno &mod)
			:
			propModule(mod),
			Cf_(mod.Cf_),
			L_(mod.L_),
			D_(mod.D_) {
	}

	Fanno::~Fanno() {
		// TODO Auto-generated destructor stub
	}

	void Fanno::changeComposition() {
		N2_->X(N1_->X());
	}

	double Fanno::calculate() {

		//Find choking conditions

		double gam = N1_->gamma();
		double M1 = N1_->M();

		double Lmax = length(M1);
		Ldiff_ = Lmax - L_;

		if (Ldiff_ <= 0) {
			throw std::range_error("Too long to not choke at exit. Reduce pipe length");
		}

		double exp1 = ((2) / (gam + 1)) * (1 + (((gam - 1) / 2) * std::pow(M1, 2)));

		N2_->X(N1_->X());

		/*N2_->p_ =  (N1_->p_ * M1) * std::pow(exp1,0.5);
        N2_->T_ =  (N1_->T_) * exp1;
        N2_->U_ =  (N1_->U_/M1) * std::pow(exp1,0.5);*/


		//double pstar = (N1_->p_ * M1) * std::pow(exp1,0.5);
		//double Tstar = (N1_->T_) * exp1;
		//double Ustar = (N1_->U_/M1) * std::pow(exp1,0.5);
		double pstar = (N1_->getPress() * M1) * std::pow(exp1, 0.5);
		double Tstar = (N1_->getTemp()) * exp1;
		double Ustar = (N1_->getU() / M1) * std::pow(exp1, 0.5);


		//Find Exit Conditions

		double Mmin, Mmax;
		if (N1_->M() <= 1) {
			Mmin = 0.0;
			Mmax = 1.0;
		} else {
			Mmin = 1.0;
			Mmax = 5.0;
		}

		double M2 = nonLinear<Fanno>::funSolve(&Fanno::lengthSolver, Mmin, Mmax, *this, 0);

		double exp2 = ((2) / (gam + 1)) * (1 + (((gam - 1) / 2) * std::pow(M2, 2)));

		//N2_->p_ =  (pstar/M2) / std::pow(exp2,0.5);
		//N2_->T_ =  Tstar/exp2;
		//N2_->U_ =  (Ustar * M2) / std::pow(exp2,0.5);
		N2_->setTP(Tstar / exp2, (pstar / M2) / std::pow(exp2, 0.5));
		N2_->setU((Ustar * M2) / std::pow(exp2, 0.5));


		return 0;

	}


	double Fanno::length(double Mach) {

		double gam = N1_->gamma();
		double M = std::pow(Mach, 2);
		double L = (D_ / (4 * Cf_)) * (((1 - M) / (gam * M)) + (((gam + 1) / (2 * gam)) * std::log(
				(M) / (((2) / (gam + 1)) * (1 + (((gam - 1) / 2) * M))))));

		return L;
	}


	double Fanno::lengthSolver(double &Mach, const void *par[]) {

		double gam = N1_->gamma();
		double M = std::pow(Mach, 2);
		double L = (D_ / (4 * Cf_)) * (((1 - M) / (gam * M)) + (((gam + 1) / (2 * gam)) * std::log(
				(M) / (((2) / (gam + 1)) * (1 + (((gam - 1) / 2) * M))))));

		return L - Ldiff_;
	}

	std::string Fanno::typeName() const {

		std::string name = "Fanno";
		return name;
	}

void Fanno::serialize(Archive& ar) const{
	propModule::serialize(ar);

	ar.put("Friction coefficient", Cf_);
	ar.put("Length", L_);
	ar.put("Hydraulic diameter", D_);
}

void Fanno::unserialize(const Archive& ar) {
	propModule::unserialize(ar);

	Cf_ = ar.getDouble("Friction coefficient");
	L_ = ar.getDouble("Length");
	D_ = ar.getDouble("Hydraulic diameter");
}
}
