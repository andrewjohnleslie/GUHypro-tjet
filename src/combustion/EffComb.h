/*!
 *  \brief      Combustion modelling with constant efficiency
 *  \details    This model assumes the combustion to be complete
 *              but it considers an efficiency applied on the reaction heat.
 *
 *              The class has a template parameter `Delta` that can be:
 *               - the class `balanceMach` if only combustion is wanted
 *               - the class `Friction` if the friction needs to be included in the model.
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

#ifndef EFFCOMB_H_
#define EFFCOMB_H_

#include "Combustion.h"

namespace hypro {
	template<class Delta>
	class EffComb : public Combustion<Delta> {
	protected:
		///Flux definition
		/**The function is reimplemented to remove the heat
         * not actually used by the combustion due to its efficiency
         * lower than one.
         */
		void delta(
				double &dG,
				double &dI,
				double &dH) const;

	public:
		///Combustion efficiency
		/**It is defined as multiplicative factor that reduces
         * the reaction heat of the complete combustion.
         */
		double eta_;

		EffComb();

		EffComb(Node *N1, Node *N2);

		EffComb(std::string name, Node *N1, Node *N2);

		EffComb(const EffComb &mod);

		virtual ~EffComb();

		propModuleSERIALtemplate(EffComb);

	virtual void serialize(Archive& ar) const;
	virtual void unserialize(const Archive& ar);


};
}
//BOOST_CLASS_EXPORT_KEY(EffComb<balanceMach>);
//BOOST_CLASS_EXPORT_KEY(EffComb<Friction>);

#endif /* EFFCOMB_H_ */
