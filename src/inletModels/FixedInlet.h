/*!
 *  \brief      Fixed geometry inlet.
 *  \details    The inlet is supposed to work supercritical,
 *              thus with a first oblique shock on the forebody and a normal shock
 *              within the inlet divergent.
 *              The position of the shock in the divergent determines the pressure drop
 *              that is chosen by means of the choking feedback coming from downstream modules.
 *              The mass flow rate is calculated knowing the coefficient of critical
 *              flow as function of free stream Mach number (see class ZitaInlet).
 *
 *              The inlet is modelled by means of a collection of modules.
 *              two internal nodes are defined:
 *              	- 0 the actual intake of the inlet duct
 *              	- 1 the throat of the inlet duct
 *              the external nodes are the free stream and the end of the inlet duct.
 *
 *              The modules used in this base class are:
 *              	- ZitaInlet between node `N1_` and internal node 0;
 *              	- EffIsenDuct between internal node 0 and `N2_`;
 *              	- node 1 is not used
 *
 *              If reduced, it decreases the pressure drop across the divergent (i.e. the EffIsenDuct module)
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

#ifndef FIXEDINLET_H_
#define FIXEDINLET_H_

#include "Inlet.h"

namespace hypro {
	class FixedInlet : public Inlet {
	public:
		FixedInlet();

		FixedInlet(Node *N1, Node *N2);

		FixedInlet(std::string, Node *N1, Node *N2);

		~FixedInlet();

	virtual void serialize(Archive& ar) const;
	virtual void unserialize(const Archive& ar);

	double calculate();

		double pdropmax_; ///<Maximum total pressure drop after the forebody

		propModuleSERIAL(FixedInlet);


		void reduce(const double &q);

		double reduce() const;

		void unreduce();

		// template<class Archive>
};
}


//BOOST_CLASS_EXPORT_KEY(FixedInlet);



#endif /* FIXEDINLET_H_ */
