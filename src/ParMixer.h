/*!
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

#ifndef PARMIXER_H_
#define PARMIXER_H_

#include "PropModule3.h"


namespace hypro {
///Model of mixing between two flows
/**The model has two inputs nodes `N1_` and `N3_`
 */
	class ParMixer : public PropModule3 {
	protected:
		///check that the Areas constrains of the module are satisfied
		/**The area of the output node is imposed equal to the sum of areas of the two input nodes.
         * In case area of node N2_ is limited, the range error of Node::A method will stop the calculation.
         * Moreover after the calculation of the module the Area of N2_ cannot be changed,
         * otherwise the module calculation would become invalid.
         * For this reason the limits of Area for N2_ are set to avoid changes.
         * This method can be reimplemented in case the area constraints are different.
         */
		void checkAreas();

		void changeComposition();

		void delta(
				double &dG,
				double &dI,
				double &dH) const;

	public:
		double teta_; ///<Angle of the flow from node `N3_` with respect to the flow of node `N1_`
		double eff_; ///<Mixing efficiency.

		ParMixer();

		ParMixer(std::string name, Node *N1, Node *N2, Node *N3, int ID = defID_);

		ParMixer(const ParMixer &mod);

		virtual ~ParMixer();

		virtual GPObject &duplicate() { return *(new ParMixer(*this)); }

		propModuleSERIAL(ParMixer);

	virtual void serialize(Archive& ar) const;
	virtual void unserialize(const Archive& ar);


};

//BOOST_CLASS_EXPORT_KEY(ParMixer);
}
#endif /* PARMIXER_H_ */
