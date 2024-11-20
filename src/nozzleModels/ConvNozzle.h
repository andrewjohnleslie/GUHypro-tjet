/*!
 *  \brief      Simple convergent nozzle.
 *  \details    
 *              The exit is assumed sonic and the exit area is calculated as consequence.
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

#ifndef CONVNOZZLE_H_
#define CONVNOZZLE_H_

#include "solvers/isentropicDuct.h"

namespace hypro {
	class ConvNozzle : public isentropicDuct {
	public:

	const Node* freeStream_; ///<Pointer to the free stream Node

	ConvNozzle();
	ConvNozzle(std::string name, Node* N1, Node* N2);
	ConvNozzle(std::string name, Node* N1, Node* N2, int ID);
	ConvNozzle(const ConvNozzle& mod);
	virtual ~ConvNozzle();

	virtual void serialize(Archive& ar) const;
	virtual void unserialize(const Archive& ar);


		propModuleSERIAL(ConvNozzle);

		double calculate();

		GPObject &duplicate() { return *(new ConvNozzle(*this)); }

		virtual bool canChoke()const;
	};
}
//BOOST_CLASS_EXPORT_KEY(hypro::ConvNozzle);
#endif /* CONVNOZZLE_H_ */
