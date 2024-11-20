/*!
 *  \brief      Model of a constant area duct with friction
 *  \details    The friction is modelled with a constant friction coefficient.
 *              The model uses the balanceMach solver so it is correct also for
 *              real gases and can be coupled with other models.
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

#ifndef FRICTION_H_
#define FRICTION_H_

#include "balanceMach.h"

namespace hypro {
	class Friction : public balanceMach {
	protected:
		Friction(std::string name, Node *in, Node *N2, int ID, int args);

		Friction(std::istream &is, const NodeMap &nodeMap, const ModuleMap &moduleMap);

		void changeComposition();

		void delta(
				double &dG,
				double &dI,
				double &dH) const;

	public:
		double Cf_; ///<Friction coefficient
		double L_; ///<Length of the duct
		double D_; ///<Hydraulic diameter of the duct
		double heat_flux; ///<Heat Flux Flag
		double enthalpy_wall;
		double stanton_number;
		double qdot;

		Friction();

		Friction(std::string name, Node *N1, Node *N2, int ID = defID_);

		Friction(const Friction &mod);

		virtual ~Friction();

		virtual void serialize(Archive& ar) const;
	virtual void unserialize(const Archive& ar);

		virtual std::string typeName() const;
	};
}

//BOOST_CLASS_EXPORT_KEY(Friction);

#endif /* FRICTION_H_ */
