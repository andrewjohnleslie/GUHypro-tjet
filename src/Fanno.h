/*!
 *  \brief      Fanno flow model
 *  \details    Valid for a constant area duct.
 *              Strictly valid only for constant gamma gases.
 *  \author     Mark De Luca
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

#ifndef FANNO_H_
#define FANNO_H_

#include "core/propModule.h"

namespace hypro {
	class Fanno : public propModule {
	protected:

		Fanno(std::string name, Node *N1, Node *N2, int ID, int args);

		///Set the composition of node `N2_`
		/**Not used anywhere.
         * @todo not sure why needed, it should be removed.
         */
		void changeComposition();

	public:
		double Cf_; ///<Friction coefficient
		double L_; ///<Length of the duct
		double D_; ///<Hydraulic diameter of the duct

		///Difference between sonic lenght and actual length
		/**@todo: Here not sure why it is needed.
         * 	However the variable should not be public.
         */
		double Ldiff_;

		Fanno(Node *N1, Node *N2);

		Fanno(std::string name, Node *N1, Node *N2);

		Fanno(std::istream &is, const NodeMap &nodeMap, const ModuleMap &moduleMap);

		Fanno(const Fanno &mod);

		virtual ~Fanno();

		double calculate();

		///<Fanno length equation
		/**It is the equation at the basis of the Fanno flow model.
         * @todo this function has no reason to be public
         *
         * @param Mach Mach number
         * @return length of the duct that would produce sonic outlet starting from Mach
         */
		double length(double Mach);

		///<Version of length adapted for non linear solver
		/**
         * @param Mach mach number
         * @param par parameter vector. Empty in this case.
         * @return residual
         */
		double lengthSolver(double &Mach, const void *par[]);

		virtual void serialize(Archive& ar) const;
		virtual void unserialize(const Archive& ar);

		virtual std::string typeName() const;
	};
}

#endif /* FANNO_H_ */
