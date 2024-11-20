/*!
 *  \brief      Model of mixing between two flows
 *  \details    The model has two inputs nodes `N1_` and `N3_`
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

#ifndef MIXER_H_
#define MIXER_H_

#include "PropModule3.h"

namespace hypro {
	class Mixer : public PropModule3 {
	protected:

		void changeComposition();

		void delta(
				double &dG,
				double &dI,
				double &dH) const;

	public:
		double teta_; ///<Angle of the flow from node `N3_` with respect to the flow of node `N1_`
		double eff_; ///<Mixing efficiency.

		Mixer();

		Mixer(std::string name, Node *N1, Node *N2, Node *N3, int ID = defID_);

		Mixer(const Mixer &mod);

		virtual ~Mixer();

		virtual GPObject &duplicate() { return *(new Mixer(*this)); }

		propModuleSERIAL(Mixer);

	virtual void serialize(Archive& ar) const;
	virtual void unserialize(const Archive& ar);

};
}

#endif /* MIXER_H_ */
