/*!
 *  \brief      Fuel injection determined by equivalent ratio
 *  \details    The amount of fuel injected is calculated based on the
 *              provided value of `PhiMax_` and the quantity of oxidizer present in the main flow.
 *              In case the module need to be reduced the actual `Phi_` is diminished.
 *              In case the fuel injected would be greater than `mfrMax_`, the actual `Phi_` is
 *              also reduced to do not overtake the threshold.
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

#ifndef INJECTIONPHINOTEMP_H_
#define INJECTIONPHINOTEMP_H_

#include <core/propModule.h>

namespace hypro {
	class balanceMach;

	class Friction;

	template<class Delta>
	class InjectionPhiNoTemp : public Delta {

	protected:
		void changeComposition();

		void delta(
				double &dG,
				double &dI,
				double &dH) const;

	public:

		///Stoichiometric coefficient of the reactants
		/**These are the coefficients used to write the left member of the reaction formula.
         * Rcoeff_[i] is the coefficient of specie i-th in the gas model referenced by the nodes.
         * Rcoeff_[i]=0 indicates that the specie i-th does not participate to the reaction.
         */
		std::vector<double> Rcoeff_;


		double Phi_;    ///<Actual Equivalent ratio
		double PhiMax_; ///<Maximum Equivalent ratio
		double mfrMax_; ///<Maximum mass flow rate of fuel
		std::vector<double> Xf_; ///<Injected fuel composition
		double Tf_; ///<Injected fuel temperature

		InjectionPhiNoTemp();

		InjectionPhiNoTemp(std::string name, Node *in, Node *N2);

		InjectionPhiNoTemp(std::string name, Node *in, Node *N2, int ID);

		InjectionPhiNoTemp(const InjectionPhiNoTemp &mod);

		virtual ~InjectionPhiNoTemp();

		virtual GPObject &duplicate() { return *(new InjectionPhiNoTemp(*this)); }

		propModuleSERIALtemplate(InjectionPhiNoTemp);

		double calculate();

		///Finds the stoichiometric composition
		/**
         * @return molar ratio between fuel and oxidizer in stochiometric conditions
         */
		double findEqv() const;

		void reduce(const double &Phi);

		double reduce() const;

		bool isreduceable() const;

		void unreduce();

		std::vector<double> propellantMfr() const;

	virtual void serialize(Archive& ar) const;
	virtual void unserialize(const Archive& ar);


};

//BOOST_CLASS_EXPORT_KEY(InjectionPhiNoTemp<balanceMach>);
//BOOST_CLASS_EXPORT_KEY(InjectionPhiNoTemp<Friction>);
}

#endif /* INJECTIONPHINOTEMP_H_ */
