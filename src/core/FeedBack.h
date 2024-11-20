/*!
 *  \brief			Class used to define chocking feedbacks in the optimizer
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

#ifndef FEEDBACK_H_
#define FEEDBACK_H_

#include "MyNode.h"
#include <core/Collection.h>
namespace hypro {
	class FeedBack : public MyNode {

	public:
		const int layer_; ///<Number defining the number of modules to count before defining a the feedback destination
		Collection::ModulePtr from_ = NULL; ///<pointer to the feedback origin
		Collection::ModulePtr to_ = NULL; ///<Pointer to the feedback destination

		FeedBack(std::string name, int layer, int ID = defID_);

		FeedBack(const FeedBack &feedB);

		~FeedBack();

		GPObject &duplicate() {
			return *(new FeedBack(*this));
		}
	};

	inline bool operator==(const FeedBack &lhs, const FeedBack &rhs) { return lhs.layer_ == rhs.layer_; }

	inline bool operator!=(const FeedBack &lhs, const FeedBack &rhs) { return !(lhs == rhs); }
}
#endif /* FEEDBACK_H_ */
