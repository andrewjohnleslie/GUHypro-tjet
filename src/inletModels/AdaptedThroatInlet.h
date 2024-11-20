/*!
 *  \brief      Inlet with variable geometry assumed always adapted, changing the size of the throat
 *  \details    The inlet is modelled by means of a collection of modules.
 *              two internal nodes are defined:
 *              	- 0 the actual intake of the inlet duct
 *              	- 1 the throat of the inlet duct
 *              the external nodes are the free stream and the end of the inlet duct.
 *
 *              The modules used in this base class are:
 *              	- chokedConv between node `N1_` and internal node 1;
 *              	- NeutralLink, substituted by an isentropicDuct if adaption limits are reached, between node `N1_` and 0;
 *              	- SupersonicInlet between internal node 0 and `N2_`.
 *
 *              If reduced, it diminishes the throat area (i.e. `nodes_[1]->A()`)
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

#ifndef ADAPTEDTHROATINLET_H_
#define ADAPTEDTHROATINLET_H_

#include "Inlet.h"

namespace hypro {
	class AdaptedThroatInlet : public Inlet {
	protected:
		using Inlet::duplicate;
	public:
		AdaptedThroatInlet();

		AdaptedThroatInlet(std::string name, Node *N1, Node *N2, int ID = defID_);

		~AdaptedThroatInlet();

		virtual GPObject &duplicate();

		propModuleSERIAL(AdaptedThroatInlet);

		double calculate();

		void reduce(const double &A);

		double reduce() const;

		bool isreduceable() const;

		void unreduce();

		virtual Glib::RefPtr<ModuleGraph> draw();
};

}
//BOOST_CLASS_EXPORT_KEY(AdaptedThroatInlet);

#endif /* ADAPTEDTHROATINLET_H_ */
