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

#include "MultiModeMachP0.h"

namespace hypro {
	MultiModeMachP0::MultiModeMachP0(std::string Name, Node *N1, Node *N2)
			:
			MultiModeMach(Name, N1, N2) {
	}

	MultiModeMachP0::~MultiModeMachP0() {
		// TODO Auto-generated destructor stub
	}

	void MultiModeMachP0::addRange() {
		std::array<double, 2> defRange;
		defRange[0] = 0.0;
		defRange[1] = std::numeric_limits<double>::infinity();
		rangesP0_.push_back(defRange);

		MultiModeMach::addRange();
	}

	int MultiModeMachP0::autoSelect() const {
		std::list<int> p = priority();
		int modeId = -1;
		for (std::list<int>::iterator it = p.begin(); it != p.end(); ++it) {
			if ((N1_->M() >= rangesM_[*it][0] && N1_->M() <= rangesM_[*it][1]) &&
				(N1_->p0() >= rangesP0_[*it][0] && N1_->p0() <= rangesP0_[*it][1])) {
				modeId = *it;
				break;
			}
		}

		if (modeId == -1) {
			std::ostringstream er;
			er << "Error: Mode not found for conditions M=" << N1_->M() << " and p0=" << N1_->p0();
			throw std::runtime_error(er.str());
		}

		return modeId;
	}

	double MultiModeMachP0::events() const {
		double p0 = N1_->p0();
		double e = std::max(rangesP0_[presentMode_][0] - p0,
							p0 - rangesP0_[presentMode_][1]) / p0;
		return std::max(MultiModeMach::events(), e);
	}
}
