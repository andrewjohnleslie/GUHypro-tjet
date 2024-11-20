/*!
 *  \brief      Variable area duct with efficiency calculated by means of MIL standard schedule with Mach number
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

#ifndef SUPERSONICINLET_H_
#define SUPERSONICINLET_H_

#include "solvers/EffIsenDuct.h"

namespace hypro {
	class SupersonicInlet : public EffIsenDuct {
		///Mach number points
		/**
         * @deprecated not used in the current version
         */
		static const std::vector<double> Mfit_;

		///Pressure drop values
		/**
         * @deprecated not used in the current version
         */
		static const std::vector<double> pfit_;

	public:
		SupersonicInlet();  ///<Empty constructor
		SupersonicInlet(std::string name, Node *in, Node *N2, int ID = defID_);

		virtual ~SupersonicInlet();

		virtual GPObject &duplicate() { return *(new SupersonicInlet(*this)); }

		propModuleSERIAL(SupersonicInlet);

		virtual void pressLoss();

		double calculate();

	virtual void serialize(Archive&)const;

//	template<class Archive>
//	void serialize(Archive & ar, const unsigned int version);
};

//BOOST_CLASS_EXPORT_KEY(SupersonicInlet);
}
#endif /* SUPERSONICINLET_H_ */
