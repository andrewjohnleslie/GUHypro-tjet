/*!
 *  \brief      Variable geometry nozzle always adapted
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

#ifndef ADAPTEDNOZZLE_H_
#define ADAPTEDNOZZLE_H_

#include <gp.h>
#include <core/Node.h>
#include "core/propModule.h"

namespace hypro {
	class isentropicDuct;

	class EffIsenDuct;

	template<class Eff>
	class AdaptedNozzle : public Eff {
	public:
		Node *freeStream_; ///<Pointer to the free stream Node

		AdaptedNozzle();

		AdaptedNozzle(std::string name, Node *N1, Node *N2);

		AdaptedNozzle(std::string name, Node *N1, Node *N2, int ID);

		AdaptedNozzle(const AdaptedNozzle &mod);

		virtual ~AdaptedNozzle();

		virtual GPObject &duplicate() { return *(new AdaptedNozzle(*this)); }

		propModuleSERIALtemplate(AdaptedNozzle);

		virtual double calculate();

	virtual void serialize(Archive& ar) const;
	virtual void unserialize(const Archive& ar);	};
}

//BOOST_CLASS_EXPORT_KEY(AdaptedNozzle<isentropicDuct>);
//BOOST_CLASS_EXPORT_KEY(AdaptedNozzle<EffIsenDuct>);

#endif /* ADAPTEDNOZZLE_H_ */
