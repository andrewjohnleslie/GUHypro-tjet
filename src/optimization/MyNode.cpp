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

#include "MyNode.h"

namespace hypro {
MyNode::MyNode(int nVal, std::string str, Type type, int args)
:
		GPNode(nVal,str.c_str(),args),
		type_(type),
		Name_(str){
}

MyNode::MyNode(const MyNode& gpo)
:
		GPNode(gpo),
		type_(gpo.type_),
		Name_(gpo.Name_){
}

	MyNode::~MyNode() {
		// TODO Auto-generated destructor stub
	}

}
