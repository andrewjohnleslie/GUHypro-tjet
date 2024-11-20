/*!
 *  \brief      Isentropic duct with efficiency dependent on inlet Mach
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

#ifndef EFFMACHISENDUCT_H_
#define EFFMACHISENDUCT_H_

#include "isentropicDuct.h"
#include "interpolation.h"

namespace hypro {
	class EffMachIsenDuct : public isentropicDuct {
	protected:
		virtual void pressLoss();

		virtual void thermLoss();

	public:
		///Total pressure efficiency as function of inlet Mach
		/**
         * etap_ is given as a set of points by the user.
         * etap_[0] is Mach and etap_[1] is the actual efficiency.
         * if etap_ is empty the efficiency is assumed 1 by default.
         */
		std::array<std::vector<double>, 2> etap_;

		///Total temperature efficiency as function of inlet Mach
		/**
         * etaT_ is given as a set of points by the user.
         * etaT_[0] is Mach and etaT_[1] is the actual efficiency.
         * if etaT_ is empty the efficiency is assumed 1 by default.
         */
		std::array<std::vector<double>, 2> etaT_;

		EffMachIsenDuct();

		EffMachIsenDuct(Node *N1, Node *N2);

		EffMachIsenDuct(std::string name, Node *N1, Node *N2);

		EffMachIsenDuct(const EffMachIsenDuct &mod); //copy constructor
		virtual ~EffMachIsenDuct();

		propModuleSERIAL(EffMachIsenDuct);

		virtual void serialize(Archive& ar) const;
		virtual void unserialize(const Archive& ar);
};
}
//BOOST_CLASS_EXPORT_KEY(hypro::EffMachIsenDuct);


#endif /* EFFMACHISENDUCT_H_ */
