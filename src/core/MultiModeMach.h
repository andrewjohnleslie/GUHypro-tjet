/*!
 *  \brief      Multi mode engine with switching based on free stream Mach number
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

#ifndef MULTIMODEMACH_H_
#define MULTIMODEMACH_H_

#include "MultiMode.h"

namespace hypro {
	class MultiModeMach : public MultiMode {
	protected:

		///Add an infinite range at the end of the vector `rangeM_`
		virtual void addRange();

	public:
		///Mach ranges of operations for each mode
		/**For each mode the range is defined by means of a low and high
         * Mach number thresholds.
         * Thus the vector has to have length equal to the number of modes.
         * The vector is update with a default infinite range whenever a new mode is added.
         */
		std::vector<std::array<double, 2> > rangesM_;

	MultiModeMach(){}

		MultiModeMach(std::string Name, Node *N1, Node *N2);

		~MultiModeMach();

		virtual GPObject &duplicate() { return *(new MultiModeMach(*this)); }

		propModuleSERIAL(MultiModeMach);

		virtual systemModule &add(std::string Name);

		virtual systemModule &add(std::string Name, int i);

		virtual int autoSelect() const;

		virtual double events() const;
	};
}
#endif /* MULTIMODEMACH_H_ */
