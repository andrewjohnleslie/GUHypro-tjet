/*!
 *  \brief      Combustion modelling with efficiency as function of phi
 *  \details    This model assumes the combustion to be complete
 *              but it considers an efficiency applied on the reaction heat.
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

#ifndef EFFPHICOMB_H_
#define EFFPHICOMB_H_

#include "combustion/Combustion.h"

namespace hypro {
	template<class Delta>
	class EffPhiComb : public Combustion<Delta> {
	protected:
		/**The function is reimplemented to remove the heat
         * not actually used by the combustion due to its efficiency
         * lower than one.
         */
		void delta(
				double &dG,
				double &dI,
				double &dH) const;

	public:

		EffPhiComb();

		EffPhiComb(Node *N1, Node *N2);

		EffPhiComb(std::string name, Node *N1, Node *N2);

		EffPhiComb(const EffPhiComb &mod);

		virtual ~EffPhiComb();

		///Efficiency as function of phi
		/**
         * eta_ is given as a set of points by the user.
         * eta_[0] is phi and zita_[1] is the actual efficiency eta.
         */
		std::array<std::vector<double>, 2> eta_;

		///Index of fuel in the gas model
		/***It is used to define the phi
         * so that the phi is the ratio between fuel and oxider and
         * not the other way around.
         */
		std::size_t iFuel_;

		propModuleSERIALtemplate(EffPhiComb);

	virtual void serialize(Archive& ar) const;
	virtual void unserialize(const Archive& ar);


};
}
//BOOST_CLASS_EXPORT_KEY(EffPhiComb<balanceMach>);
//BOOST_CLASS_EXPORT_KEY(EffPhiComb<Friction>);
#endif /* EFFCOMB_H_ */
