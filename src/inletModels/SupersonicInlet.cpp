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

#include "SupersonicInlet.h"
#include "interpolation.h"

namespace hypro {
	const std::vector<double> SupersonicInlet::Mfit_{5.0, 8.0, 9.0, 10.0};
	const std::vector<double> SupersonicInlet::pfit_{0.71, 0.32, 0.21, 0.13};

	SupersonicInlet::SupersonicInlet(std::string name, Node *N1, Node *N2, int ID)
			:
			EffIsenDuct(name, N1, N2, ID) {

	}

	SupersonicInlet::SupersonicInlet() :
			EffIsenDuct() {
	}

	SupersonicInlet::SupersonicInlet(std::istream &is, const NodeMap &nodeMap, const ModuleMap &moduleMap)
			:
			EffIsenDuct(is, nodeMap, moduleMap) {
	}

	SupersonicInlet::~SupersonicInlet() {
		// TODO Auto-generated destructor stub
	}

	void SupersonicInlet::pressLoss() {
		double M1 = N1_->M();
		if (M1 <= 1.0) {
			//The constant pressure efficiency is used only in subsonic
			EffIsenDuct::pressLoss();
		} else if (M1 > 1 && M1 < 7) {
			//p02_p01_ = 1.0 - 0.0766*std::pow(M1 - 1.0, 1.35);
//			p02_p01_ = 1.0 - 0.075 * std::pow(M1 - 1.0, 1.35);
			p02_p01_ = 0.937;
		} else if (M1 >= 7) {
			p02_p01_ = 0.341;
		}
	}

// Reimplements calculate() method to make sure choked_ is set properly
	double SupersonicInlet::calculate() {
		choked_ = N1_->M() >= 1.0;
		return EffIsenDuct::calculate();
	}


void SupersonicInlet::serialize(Archive& ar)const{
	EffIsenDuct::serialize(ar);
}

//BOOST_CLASS_EXPORT_IMPLEMENT(SupersonicInlet);
}
