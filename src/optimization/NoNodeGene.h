/*!
 *  \brief      GP gene to define a scalar parameter
 *  \author     Alessandro Mogavero
 *  \n
 *  \email      alessandro.mogavero@strath.ac.uk
 *  \version    1.0
 *  \copyright  Copyright 2017 Alessandro Mogavero
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

#ifndef NONODEGENE_H_
#define NONODEGENE_H_

#include <MyGene.h>
#include <systemModule.h>

namespace hypro{
class GPTypedGeneSet;

class NoNodeGene : public MyGene
{
	int value_;

public:

	NoNodeGene():
		MyGene(),
		value_(0){}

	///Constructor based on GP node
	NoNodeGene (int value);

	///Copy constructor
	NoNodeGene (const NoNodeGene& gpo);

	virtual ~NoNodeGene();

	virtual MyNode::Type nodeType()const;

	virtual propModule& getNode()const;

	virtual int value()const;

	char* load (istream& is);

	virtual MyGene& operator= (const MyGene& gpo);

	virtual std::string name()const;
};
}
#endif /* NONODEGENE_H_ */
