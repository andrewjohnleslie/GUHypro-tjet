/*!
 *  \author     Robert Garner
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

#include "ScramjetInlet.h"
#include "interpolation.h"

namespace hypro {
    const std::vector<double> ScramjetInlet::Mfit_{5.0, 8.0, 9.0, 10.0};
    const std::vector<double> ScramjetInlet::pfit_{0.71, 0.32, 0.21, 0.13};

    ScramjetInlet::ScramjetInlet():
        isentropicDuct(){
        }

    ScramjetInlet::ScramjetInlet(Node *N1, Node *N2)
            :
            isentropicDuct("inlet", N1, N2) {

    }

    ScramjetInlet::ScramjetInlet(std::string name, Node *N1, Node *N2)
            :
            isentropicDuct(name, N1, N2) {

    }

    ScramjetInlet::ScramjetInlet(std::istream &is, const NodeMap &nodeMap, const ModuleMap &moduleMap)
            :
            isentropicDuct(is, nodeMap, moduleMap) {
    }

    ScramjetInlet::~ScramjetInlet() {
        // TODO Auto-generated destructor stub
    }

    void ScramjetInlet::pressLoss() {
        double M1 = N1_->M();
        if (M1 < 5.0) {
            p02_p01_ = 1;
        } else if (M1 >= 5) {
            p02_p01_ = 0.00232359 * std::pow(M1, 3) - 0.0479577 * std::pow(M1, 2) + 0.19649942 * M1 + 0.64185476;
        }
    }

// Reimplementats calculate() method to make sure choked_ is set properly
    double ScramjetInlet::calculate() {
        choked_ = N1_->M() >= 1.0;
        choked_ = false;
        return isentropicDuct::calculate();
    }

    double ScramjetInlet::drag() const {
        //here we should add also the non adaption drag
        //This definition is good for optimization but not for normal calculations
        return N1_->A() * N1_->rho() * std::pow(N1_->getU(), 2.0);
    }

    void ScramjetInlet::serialize(Archive& ar)const{
	isentropicDuct::serialize(ar);
}

}