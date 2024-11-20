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

#ifndef CHOKEDCONV_H_
#define CHOKEDCONV_H_

#include "SupersonicInlet.h"

namespace hypro {
	class chokedConv :
			public SupersonicInlet {
	public:
		enum ThroatArea {
			ASSIGNED, CALCULATED
		};
		///Changes the behaviour of the module
		/**ASSIGNED means the area of node `N1_` will be adapted to ensure sonic conditions
         * CALCULATED means the `N2_` area will be adapted to ensure sonic conditions
         */
		ThroatArea throat_;

		/// Value of free stream Mach used in case of engine fixed point condition
		static double startToll_;

		chokedConv();  ///<Empty constructor
		chokedConv(std::string name, Node *N1, Node *N2, ThroatArea throat = CALCULATED, int ID = defID_);

		chokedConv(const chokedConv &mod);

		virtual ~chokedConv();

		virtual GPObject &duplicate() { return *(new chokedConv(*this)); }

		propModuleSERIAL(chokedConv);

	virtual void serialize(Archive& ar) const;
	virtual void unserialize(const Archive& ar);

	virtual double calculate();


	virtual bool canChoke()const;


};
}
//BOOST_CLASS_EXPORT_KEY(chokedConv);


#endif /* INLET_H_ */
