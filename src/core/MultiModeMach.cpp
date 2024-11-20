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

#include "MultiModeMach.h"

namespace hypro {
	MultiModeMach::MultiModeMach(std::string Name, Node *N1, Node *N2)
			:
			MultiMode(Name, N1, N2) {
	}

	MultiModeMach::~MultiModeMach() {
		// TODO Auto-generated destructor stub
	}

	systemModule &MultiModeMach::add(std::string Name) {
		addRange();

		return MultiMode::add(Name);
	}

	systemModule &MultiModeMach::add(std::string Name, int i) {
		addRange();

		return MultiMode::add(Name, i);
	}

	void MultiModeMach::addRange() {
		std::array<double, 2> defRange;
		defRange[0] = 0.0;
		defRange[1] = std::numeric_limits<double>::infinity();
		rangesM_.push_back(defRange);
	}

	int MultiModeMach::autoSelect() const {
		std::list<int> p = priority();
		int modeId = -1;
		for (std::list<int>::iterator it = p.begin(); it != p.end(); ++it) {
			if (N1_->M() >= rangesM_[*it][0] && N1_->M() <= rangesM_[*it][1]) {
				modeId = *it;
				break;
			}
		}

		return modeId;
	}

	double MultiModeMach::events() const {
		return std::max(rangesM_[presentMode_][0] - N1_->M(),
						N1_->M() - rangesM_[presentMode_][1]);
	}
}