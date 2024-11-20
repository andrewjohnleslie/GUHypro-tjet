/*!
 *  \brief      Provides the solver for the balance equations
 *  \details    This class solves the balance equations.
 *              It has to be inherited by those classes that need to solve the balance equations
 *              of a duct with constant area.
 *              To define the problem the methods changeComposition and delta
 *              need to be reimplemented.
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

#ifndef BALANCEMACH_H_
#define BALANCEMACH_H_

#include "core/propModule.h"

//#define BOOST_MATH_INSTRUMENT
#include <boost/math/tools/roots.hpp>

namespace hypro {
	class balanceMach : public propModule {

		///Balance equation to be enforced
		/**Calculate in sequence the 3 balance equations.
         *
         * @param[in] M2 Mach number of node `N2_`
         * @param[in] par array of parameters. Empty in this case
         * @return the residual
         */
		double balanceG(double &M2, const void *par[]);

		///Energy equation
		/**Used to calculate the static temperature knowing total enthalpy and
         * Mach number.
         *
         * @param[in] T static temperature.
         * @param[in] par array of parameters.
         *  - par[0] Node for which the calculation is performed
         *  - par[1] Mach number
         *  - par[2] Total enthalpy
         */
		double THequation(double &T, const void *par[]);

		///Returns the choking residual
		/**This method is called at the beginning of the calculate() method
         * in order to check if the module is overchoked or not.
         *
         * @return the residual of the sonic case. Negative if over choked.
         */
		double choked();

	protected:
		double G1_; ///<Mass flow of node `N1_`
		double I1_; ///<Impulse of node `N1_`
		double H1_; ///<Enthalpy flow of node `N1_`

		///Assign a composition to output node `N2_`
		/**This method needs to be reimplemented in subclasses.
         * It is called before the solution of balance equation problem.
         */
		virtual void changeComposition()=0;

		///Defines fluxes
		/**This method needs to be reimplemented in subclasses.
         * The default implementation defines null fluxes.
         * It is called during the iterations of the solution of balance equation problem.
         *
         * @param[out] dG is the flux of mass
         * @param[out] dI is the flux of momentum
         * @param[out] dH is the flux of energy
         */
		virtual void delta(
				double &dG,
				double &dI,
				double &dH) const;

		///check that the Areas constrains of the module are satisfied
		/**The area of the output node is imposed equal to the area of input node.
         * In case area of node N2_ is limited, the range error of Node::A method will stop the calculation.
         * Moreover after the calculation of the module the Area of N2_ cannot be changed,
         * otherwise the module calculation would become invalid.
         * For this reason the limits of Area for N2_ are set to avoid changes.
         * This method can be reimplemented in case the area constraints are different.
         */
		virtual void checkAreas();

		balanceMach(std::string name, Node *N1, Node *N2, int ID, int args);

		balanceMach(std::istream &, const NodeMap &nodeMap, const ModuleMap &moduleMap);

	public:
		balanceMach();

		balanceMach(std::string name, Node *N1, Node *N2, int ID = defID_);

		balanceMach(const balanceMach &name);

		virtual ~balanceMach();

		virtual double calculate();

		bool solver_verbosity = false;


		using propModule::serialize;
	};
}

//BOOST_CLASS_EXPORT_KEY(hypro::balanceMach);

#endif /* BALANCEMACH_H_ */
