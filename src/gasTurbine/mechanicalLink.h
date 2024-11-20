/*!
// Created by Mark De Luca on 07/12/18.
//
// For any queries, contact m.de-luca.1@research.gla.ac.uk
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
 * along with HyPro.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef MECHANICALLINK_H_
#define MECHANICALLINK_H_

#include <string>

namespace hypro {

	class mechanicalLink {

		static std::string defName_;
		static unsigned nextID_;
		double Power_{0};


	public:

		unsigned int ID_;

		double Cp_, Cp2_, kH_, mfr_, Ht_;
		double eff_mech_;

		std::string Name_;

		mechanicalLink();
		mechanicalLink(const mechanicalLink& Link);
		~mechanicalLink();

		double getPower();
		void setPower(double P);
	};
}

#endif
