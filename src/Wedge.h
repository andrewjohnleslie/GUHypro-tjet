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

#ifndef WEDGE_H_
#define WEDGE_H_

#include "core/propModule.h"

namespace hypro {
	class Wedge : public propModule {
		///Oblique shock equation
		/**Evaluates the oblique shock equations in case of constant gamma gas.
         *
         * @param[in] M1 Mach number before the shock
         * @param[in] gamma specific heats ratio
         * @param[in] delta wedge angle [deg]
         * @param[out] M2 Mach number after the shock
         * @param[out] p2_p1 pressure ratio before and after the shock
         * @param[out] T2_T1 temperature ratio before and after the shock
         */
		void obliqueShock(
				const double &M1, const double &gamma, const double &delta,
				double &M2, double &p2_p1, double &T2_T1);

		///Equation to be solved for the oblique shock problem
		/**
         * @param eps Angle of the shock with respect the stream before the shock
         * @param par parameter vector
         * 	- par[0] Mach before the shock
         * 	- par[1] gamma, specific heats ratio
         * 	- par[2] delta, wedge angle
         * @return residual
         */
		double deltaBeta(double &eps, const void *par[]);

		///Normal shock equation
		/**Evaluates the normal shock equations in case of constant gamma gas.
         *
         * @param[in] M1 Mach number before the shock
         * @param[in] gamma specific heats ratio
         * @param[out] M2 Mach number after the shock
         * @param[out] p2_p1 pressure ratio before and after the shock
         * @param[out] T2_T1 temperature ratio before and after the shock
         */
		static void shock(
				const double &M1, const double &gamma,
				double &M2, double &p2_p1, double &T2_T1);

	public:
		double delta_; ///<Wedge angle `[deg]`

	Wedge(){}

		Wedge(std::string name, Node *N1, Node *N2, int ID = defID_);

		Wedge(const Wedge &mod);

		~Wedge();
	virtual void serialize(Archive& ar) const;
	virtual void unserialize(const Archive& ar);

	virtual double calculate();

	virtual bool canChoke()const;
		virtual GPObject &duplicate() { return *(new Wedge(*this)); }
		propModuleSERIAL(Wedge);
	};
}
#endif /* WEDGE_H_ */
