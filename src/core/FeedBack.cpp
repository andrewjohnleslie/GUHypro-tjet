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

#include "FeedBack.h"

namespace hypro {
	FeedBack::FeedBack(std::string name, int layer, int ID)
			:
			MyNode(ID, name, FEEDBACK, 0),
			layer_(layer) {
	}

	FeedBack::FeedBack(const FeedBack &feedB)
			:
			MyNode(feedB),
			from_(feedB.from_),
			layer_(feedB.layer_),
			to_(feedB.to_) {
	}

	FeedBack::~FeedBack() {
		// TODO Auto-generated destructor stub
	}
}