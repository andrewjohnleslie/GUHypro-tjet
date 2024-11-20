/*!
 *  \brief      Inlet forebody defined by means of coefficient of critical flow zita.
 *  \details    The inlet is supposed to work supercritical.
 *              The mass flow rate is calculated knowing the coefficient of critical
 *              flow as function of free stream Mach number.
 *              The maximum area of node N1 must be set equal to the maximum capture area,
 *              in a center body inlet it is the total cross section area in correspondence of the
 *              cowl (not the actual inlet area, but inlet area + cross section of center body).
 *              The maximum capture area is also the capture area in adapted conditions, so
 *              when the front shock impinges the cowl.
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

#ifndef ZITAINLET_H_
#define ZITAINLET_H_

#include "solvers/EffMachIsenDuct.h"

namespace hypro {
	class ZitaInlet : public EffMachIsenDuct {
	public:
		ZitaInlet();

		ZitaInlet(Node *N1, Node *N2);

		ZitaInlet(std::string, Node *N1, Node *N2);

		~ZitaInlet();

		///<Coefficient of critical flow as function of Mach
		/**Ratio between the actual mfr passing through the inlet
         * and the mfr that would pass if the inlet was adapted
         * (i.e. the front shock impinges the cowl).
         * zita_ is given as a set of points by the user.
         * zita_[0] is the Mach number and zita_[1] is the actual zita number.
         */
		std::array<std::vector<double>, 2> zita_;

		propModuleSERIAL(ZitaInlet);

	virtual void serialize(Archive& ar) const;
	virtual void unserialize(const Archive& ar);

	double calculate();


};
}
//BOOST_CLASS_EXPORT_KEY(ZitaInlet);


#endif /* ZITAINLET_H_ */
