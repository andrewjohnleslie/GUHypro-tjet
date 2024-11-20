/*!
 *  \brief      Variable area duct with constant temperature and pressure losses
 *  \details    This class adds constant pressure losses and temperature losses
 *              ratios to the `isentropicDuct` model.
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

#ifndef EFFISENDUCT_H_
#define EFFISENDUCT_H_

#include "isentropicDuct.h"

namespace hypro {

	class EffIsenDuct : public isentropicDuct {
	protected:
		virtual void pressLoss();

		virtual void thermLoss();

	public:
		double etap_; ///<Stagnation pressure loss
		double etaT_; ///<Stagnation temperature loss

		EffIsenDuct();  ///<Empty constructor
		EffIsenDuct(std::string name, Node *N1, Node *N2, int ID = defID_);

		EffIsenDuct(const EffIsenDuct &mod); //copy constructor
		virtual ~EffIsenDuct();

		propModuleSERIAL(EffIsenDuct);

		virtual GPObject &duplicate() { return *(new EffIsenDuct(*this)); }

		virtual void serialize(Archive& ar) const;
	virtual void unserialize(const Archive& ar);	
	};
}
//BOOST_CLASS_EXPORT_KEY(EffIsenDuct);

#endif /* EFFISENDUCT_H_ */
