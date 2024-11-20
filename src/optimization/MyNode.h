/*!
 *  \brief      GP node implementation for HyPro
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

#ifndef MYNODE_H_
#define MYNODE_H_

#include <gp.h>
#include <vector>

namespace hypro{
class MyNode: public GPNode {
protected:
	static const int defID_=0; ///<Default ID of the node. Used in constructors if omitted by the user.

public:
	///Type of the node
	const enum Type {
		MODULE = 0, //!< MODULE
		FEEDBACK = 1//!< FEEDBACK
	} type_;

	std::string Name_; ///<Name of the node/HyPro module

	MyNode():
		type_(MODULE){}

	///Constructor
	/**
	 * @param nVal node value, it identify the node during the optimization
	 * @param str  node Name
	 * @param type type of the node
	 * @param args number of node arguments (i.e. corresponding gene child).
	 */
	MyNode(int nVal, std::string str, Type type, int args=0);

	///Copy constructor
	MyNode(const MyNode& gpo);

	~MyNode();
};
}
#endif /* MYNODE_H_ */
