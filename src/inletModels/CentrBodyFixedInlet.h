/*!
 *  \brief      Inlet with central body and fixed geometry
 *  \details    The inlet is modelled by means of a collection of modules.
 *              two internal nodes are defined:
 *              	- 0 the actual intake of the inlet duct
 *              	- 1 the throat of the inlet duct
 *              the external nodes are the free stream and the end of the inlet duct.
 *
 *              The modules used in this base class are:
 *              	- Wedge between node `N1_` and internal node 0;
 *              	- EffIsenDuct between internal node 0 and `N2_`;
 *              	- internal node 1 is not used
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

#ifndef CENTRBODYFIXEDINLET_H_
#define CENTRBODYFIXEDINLET_H_

#include "Inlet.h"

namespace hypro {
	class CentrBodyFixedInlet : public Inlet {
	public:
		CentrBodyFixedInlet();

		CentrBodyFixedInlet(Node *N1, Node *N2);

		CentrBodyFixedInlet(std::string, Node *N1, Node *N2);

		~CentrBodyFixedInlet();
	virtual void serialize(Archive& ar) const;
	virtual void unserialize(const Archive& ar);


		double pdropmax_; ///<Maximum value for the pressure drop in the divergent

		propModuleSERIAL(CentrBodyFixedInlet);

		double calculate();

		void reduce(const double &q);

		double reduce() const;

		void unreduce();

		double drag() const;

		template<class Archive>
		void serialize(Archive &ar, const unsigned int version);
	};


//BOOST_CLASS_EXPORT_KEY(CentrBodyFixedInlet);
}

#endif /* FIXEDINLET_H_ */
