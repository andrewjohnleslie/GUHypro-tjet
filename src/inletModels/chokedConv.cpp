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
#include "chokedConv.h"

namespace hypro {
	double chokedConv::startToll_ = 0.001;

	chokedConv::chokedConv(std::string name, Node *N1, Node *N2, ThroatArea throat, int ID)
			:
			SupersonicInlet(name, N1, N2, ID),
			throat_(throat) {
	}

	chokedConv::chokedConv() :
			SupersonicInlet(),
			throat_(ASSIGNED) {
	}

	chokedConv::chokedConv(std::istream &is, const NodeMap &nodeMap, const ModuleMap &moduleMap)
			:
			SupersonicInlet(is, nodeMap, moduleMap) {
		std::string throat;
		is >> throat;
		if (std::strcmp(throat.c_str(), "ASSIGNED") == 0) {
			throat_ = ASSIGNED;
		} else if (std::strcmp(throat.c_str(), "CALCULATED") == 0) {
			throat_ = CALCULATED;
		} else {
			throw std::runtime_error("Error: value " + throat + " not valid for throat_ behaviour.");
		}
	}

	chokedConv::chokedConv(const chokedConv &mod)
			:
			SupersonicInlet(mod),
			throat_(mod.throat_) {
	}

	chokedConv::~chokedConv() {
		// TODO Auto-generated destructor stub
	}

	double chokedConv::calculate() {
		if (N1_->getU() == 0.0) {
			N1_->setU(startToll_ * N1_->geta());
		}

		//cout << N1_->report() << endl;
		pressLoss();
		thermLoss();
		applyLoss();
		//cout << N1_->report() << endl;
		//cout << N1mod_->report() << endl;

		double M2 = 1;
		N2_->X(N1_->X());

		double pstr_p0, Tstr_T0;
		isentropic(M2, N1_->gamma(), pstr_p0, Tstr_T0);
		double p1_p0, T1_T0;
		isentropic(N1_->M(), N1_->gamma(), p1_p0, T1_T0);

		//N2_->p_ = pstr_p0*p02_p01_*N1_->p_/p1_p0;
		//N2_->T_ = Tstr_T0*T02_T01_*N1_->T_/T1_T0;
		N2_->setTP(Tstr_T0 * T02_T01_ * N1_->getTemp() / T1_T0, pstr_p0 * p02_p01_ * N1_->getPress() / p1_p0);
		//N2_->U_ = M2*N2_->geta();
		N2_->setU(M2 * N2_->geta());
		double A2_A1 = mfrEq(M2);
		switch (throat_) {
			case ASSIGNED:
				N1_->A(N2_->A() / A2_A1);
				break;
			case CALCULATED:
				N2_->A(N1_->A() * A2_A1);
				break;
		}

		return 1.0; //never choked
	}

void chokedConv::serialize(Archive& ar) const{
	SupersonicInlet::serialize(ar);

	switch(throat_){
	case ASSIGNED:
		ar.put("Throat area calculation type", "ASSIGNED");
		break;
	case CALCULATED:
		ar.put("Throat area calculation type", "CALCULATED");
		break;
	}
}

void chokedConv::unserialize(const Archive& ar) {
	SupersonicInlet::unserialize(ar);

	std::map<std::string, ThroatArea> en_map;
	en_map["ASSIGNED"] = ASSIGNED;
	en_map["CALCULATED"] = CALCULATED;

	throat_ = en_map[ar.get<std::string>("Throat area calculation type")];
}

bool chokedConv::canChoke()const{
	return false;
}
}


//BOOST_CLASS_EXPORT_IMPLEMENT(chokedConv);
