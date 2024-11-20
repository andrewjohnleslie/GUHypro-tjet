/*!
 *  \brief      Combustion modelling
 *  \details    This model assumes the combustion to be complete.
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

#ifndef COMBUSTION_H_
#define COMBUSTION_H_

#include <solvers/balanceMach.h>

namespace hypro {
	class Friction;

	template<class Delta>
	class Combustion : public Delta {
	protected:
		void changeComposition();

	public:
		///Stoichiometric coefficient of the reactants
		/**These are the coefficients used to write the left member of the reaction formula.
         * `Rcoeff_[i]` is the coefficient of specie i-th in the gas model referenced by the nodes.
         * `Rcoeff_[i]=0` indicates that the specie i-th does not participate to the reaction.
         */
		std::vector<double> Rcoeff_;
		int Xswitch_ {0};
		std::vector<double> XSET_;

		///Stoichiometric coefficient of the products
		/**These are the coefficients used to write the right member of the reaction formula.
         * `Pcoeff_[i]` is the coefficient of specie i-th in the gas model referenced by the nodes.
         * `Pcoeff_[i]=0` indicates that the specie i-th does not participate to the reaction.
         */
		std::vector<double> Pcoeff_;

		Combustion();

		Combustion(std::string name, Node *N1, Node *N2, int ID);

		Combustion(const Combustion &mod);

		virtual GPObject &duplicate() { return *(new Combustion(*this)); }

		virtual ~Combustion();

		propModuleSERIALtemplate(Combustion);

	virtual void serialize(Archive& ar) const;
	virtual void unserialize(const Archive& ar);


};
}
//BOOST_CLASS_EXPORT_KEY(Combustion<balanceMach>);
//BOOST_CLASS_EXPORT_KEY(Combustion<Friction>);

#endif /* COMBUSTION_H_ */
