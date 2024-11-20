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
#include "ConvNozzle.h"
namespace hypro {
	ConvNozzle::ConvNozzle(std::string name, Node *N1, Node *N2)
			:
			ConvNozzle(name, N1, N2, defID_) {
	}

ConvNozzle::ConvNozzle(std::string name, Node* N1, Node* N2, int ID)
:
	freeStream_(NULL),
	isentropicDuct(name,N1,N2,ID){
}

ConvNozzle::ConvNozzle():
	freeStream_(NULL),
	isentropicDuct(){
}

ConvNozzle::ConvNozzle(std::istream& is, const NodeMap& nodeMap, const ModuleMap& moduleMap)
:
		freeStream_(NULL),
		isentropicDuct(is,nodeMap,moduleMap){
}

ConvNozzle::ConvNozzle(const ConvNozzle& mod)
:
	freeStream_(mod.freeStream_),
	isentropicDuct(mod){

	}

	ConvNozzle::~ConvNozzle() {
		// TODO Auto-generated destructor stub
	}

	double ConvNozzle::calculate() {
		double p1_p01, T1_T01;
		isentropic(N1_->M(), N1_->gamma(), p1_p01, T1_T01);

		pressLoss();
		thermLoss();

	double M2 = 1;
	applyLoss();
	N2_->A(N1_->A()*mfrEq(M2));

	if(N2_->getPress() < freeStream_->getPress()) throw std::runtime_error("Error: Convergent nozzle cannot be over-expanded.");

		return 1; //not choked
	}

bool ConvNozzle::canChoke()const{
	return false;
}

void ConvNozzle::serialize(Archive& ar) const{
	isentropicDuct::serialize(ar);

	ar.put("Free Stream node", freeStream_->ID_);
}

void ConvNozzle::unserialize(const Archive& ar) {
	isentropicDuct::unserialize(ar);

	freeStream_ = ar.getNodeRef("Free Stream node");
}
}