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

#include "interpolation.h"
#include <stdexcept>

namespace hypro {
	interpolation::interpolation() {
		// TODO Auto-generated constructor stub

	}

	interpolation::~interpolation() {
		// TODO Auto-generated destructor stub
	}

	double interpolation::linear(const std::vector<double> &X, const std::vector<double> &Y, const double &x0) {
		if (x0 < X[0] || x0 > X.back()) {
			throw std::runtime_error("Error: x0 is out of range");
		}

		double y0 = 0;
		for (std::size_t i = 1; i < X.size(); i++) {
			if (x0 <= X[i]) {
				y0 = Y[i - 1] + (x0 - X[i - 1]) * (Y[i] - Y[i - 1]) / (X[i] - X[i - 1]);
				break;
			}
		}
		return y0;
	}
}