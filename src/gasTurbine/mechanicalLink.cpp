/*!
// Created by Mark De Luca on 07/12/18.
//
// For any queries, contact m.de-luca.1@research.gla.ac.uk
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
 * along with HyPro.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "mechanicalLink.h"
#include "propModule.h"

namespace hypro {

	std::string mechanicalLink::defName_="mechanicalLink";
	unsigned mechanicalLink::nextID_ = 1;


	mechanicalLink::mechanicalLink():

			Cp_{0},
			Cp2_{0},
			Ht_{0},
			kH_{0},
			mfr_{0},
			eff_mech_{0},		
			Name_(defName_),
			ID_(nextID_){
		nextID_++;

	}

	mechanicalLink::mechanicalLink(const mechanicalLink& Link):

			Power_(Link.Power_),
			Cp_(Link.Cp_),
			Cp2_(Link.Cp2_),
			Ht_{Link.Ht_},
			kH_(Link.kH_),
			mfr_(Link.mfr_),
			eff_mech_(Link.eff_mech_),
			Name_(Link.Name_),
			ID_(nextID_){
		nextID_++;

	}

	mechanicalLink::~mechanicalLink(){
	}

	double mechanicalLink::getPower(){
		eff_mech_ = 0.95;
		double effPower_ = Power_ * (1/eff_mech_);	
		
		return effPower_;
	}

	void mechanicalLink::setPower(double P){
		Power_ += P;	
	}
}
