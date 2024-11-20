/*!
 *  \brief      Dummy module that copies `N1_` to `N2_`
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

#ifndef NEUTRALLINK_H_
#define NEUTRALLINK_H_

#include "core/propModule.h"

namespace hypro {
	class NeutralLink :
			public propModule {
	public:
		NeutralLink();

		NeutralLink(std::string name, Node *N1, Node *N2, int ID = defID_);

		NeutralLink(const NeutralLink &mod);

		~NeutralLink();

		virtual GPObject &duplicate() { return *(new NeutralLink(*this)); }

		propModuleSERIAL(NeutralLink);

		virtual double calculate();

		// template<class Archive>

		virtual bool canChoke()const;


};
}
//BOOST_CLASS_EXPORT_KEY(NeutralLink);

#endif /* NEUTRALLINK_H_ */
