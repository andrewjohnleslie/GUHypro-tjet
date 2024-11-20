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
#include "Wedge.h"
#include "nonLinear.h"

namespace hypro {
	Wedge::Wedge(std::string name, Node *N1, Node *N2, int ID)
			:
			propModule(name, N1, N2, ID),
			delta_(0.0) {
	}

	Wedge::Wedge(std::istream &is, const NodeMap &nodeMap, const ModuleMap &moduleMap)
			:
			propModule(is, nodeMap, moduleMap) {
		///serialization of delta_ is missing!
	}

	Wedge::Wedge(const Wedge &mod)
			:
			propModule(mod),
			delta_(mod.delta_) {
	}

	Wedge::~Wedge() {
		// TODO Auto-generated destructor stub
	}

	double Wedge::calculate() {
		N2_->X(N1_->X());
		if (N1_->M() > 1.0) {
			double M2, p2_p1, T2_T1;
			obliqueShock(N1_->M(), N1_->gamma(), delta_, M2, p2_p1, T2_T1);

			//N2_->p_ = p2_p1*N1_->p_;
			//N2_->T_ = T2_T1*N1_->T_;
			//N2_->U_ = M2*N2_->geta();
			N2_->setTP(T2_T1 * N1_->getTemp(), p2_p1 * N1_->getPress());
			N2_->setU(M2 * N2_->geta());
		} else {
			//N2_->p_ = N1_->p_;
			//N2_->T_ = N1_->T_;
			//N2_->U_ = N1_->U_;
			N2_->setTP(N1_->getTemp(), N1_->getPress());
			N2_->setU(N1_->getU());
		}

		return 1; //never choked
	}

	void Wedge::obliqueShock(
			const double &M1, const double &gamma, const double &delta,
			double &M2, double &p2_p1, double &T2_T1) {
		double Delta = delta * pi_ / 180;

		const void *par[3];
		par[0] = &M1;
		par[1] = &gamma;
		par[2] = &Delta;
		double r1 = deltaBeta(Delta, par);
		double r2 = r1;
		double Delta2 = Delta;
		while (r1 * r2 > 0) {
			Delta2 = Delta2 * (1.0 + 0.1);
			r2 = deltaBeta(Delta2, par);
		}
		double eps = nonLinear<Wedge>::funSolve(&Wedge::deltaBeta, Delta, Delta2, *this, par);

		double M1n = M1 * std::sin(eps);
		double M2n;
		shock(M1n, gamma, M2n, p2_p1, T2_T1);

		M2 = M2n / std::sin(eps - Delta);
	}

	double Wedge::deltaBeta(double &eps, const void *par[]) {
		const double &M1 = *(double *) par[0];
		const double &gamma = *(double *) par[1];
		const double &delta = *(double *) par[2];

		if (eps <= 0) {
			eps = 1e-10;  // avoid finding solution for eps=0
		}
		return 1 /
			   (std::tan(eps) * ((gamma + 1.0) / (2.0 * (std::pow(std::sin(eps), 2.0) - std::pow(M1, -2.0))) - 1.0)) -
			   std::tan(delta);
	}

void Wedge::shock(
		const double& M1, const double& gamma,
		double& M2, double& p2_p1, double& T2_T1){
	M2 = std::sqrt(((gamma - 1)*std::pow(M1,2.0) + 2.0)/(2.0*gamma*std::pow(M1,2.0) - (gamma - 1.0)));
	p2_p1 = (2.0*gamma*std::pow(M1,2.0) - gamma + 1.0)/(gamma + 1.0);
	T2_T1 = (2.0*gamma*std::pow(M1,2.0) - gamma + 1.0)*(gamma - 1.0 + 2.0/std::pow(M1,2.0))/std::pow(gamma + 1,2.0);
}

bool Wedge::canChoke()const{
	return false;
}

void Wedge::serialize(Archive& ar) const{
	propModule::serialize(ar);

	ar.put("Wedge angle", delta_);
}

void Wedge::unserialize(const Archive& ar) {
	propModule::unserialize(ar);

	delta_ = ar.getDouble("Wedge angle");
}
}