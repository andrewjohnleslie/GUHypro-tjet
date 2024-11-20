/*!
 *  \brief      Multi mode engine with switching based on free stream Mach number and stagnation pressure
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

#ifndef MULTIMODEMACHP0_H_
#define MULTIMODEMACHP0_H_

#include "MultiModeMach.h"

namespace hypro {
	class MultiModeMachP0 : public MultiModeMach {
	protected:
		void addRange();

	public:
		///Stagnation pressure ranges of operations for each mode
		/**For each mode the range is defined by means of a low and high
         * stagnation pressure thresholds.
         * Thus the vector has to have length equal to the number of modes.
         * The vector is update with a default infinite range whenever a new mode is added.
         */
		std::vector<std::array<double, 2>> rangesP0_;

	MultiModeMachP0(){}
	MultiModeMachP0(std::string Name, Node* N1, Node* N2);
	virtual ~MultiModeMachP0();

		virtual GPObject &duplicate() { return *(new MultiModeMachP0(*this)); }

		propModuleSERIAL(MultiModeMachP0);

		int autoSelect() const;

		virtual double events() const;
	};
}

#endif /* MULTIMODEMACHP0_H_ */
