/*!
 *  \brief      Provides functions to perform interpolations
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

#ifndef INTERPOLATION_H_
#define INTERPOLATION_H_

#include <vector>
//#include "error.H"

namespace hypro {
	class interpolation {
	public:
		interpolation();

		virtual ~interpolation();

		///Linear interpolation
		/**
         * @param X independent variable data points
         * @param Y dependent variable data points
         * @param x0 point where the interpolation is wanted
         * @return y0 the interpolated value
         */
		static double linear(const std::vector<double> &X, const std::vector<double> &Y, const double &x0);
	};
}
#endif /* INTERPOLATION_H_ */
