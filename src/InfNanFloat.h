/*
 *  \brief      Class used to serializa float number that can be inf or nan
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

#ifndef SRC_INFNANFLOAT_H_
#define SRC_INFNANFLOAT_H_

//#include <boost/serialization/split_member.hpp>

class InfNanFloat {
	double value_; ///<Actual value to be serialized
	bool isInf_;  ///<True if value is infinity
	bool isNan_;  ///<True if value is not a number
public:
	InfNanFloat(double value);
	virtual ~InfNanFloat();

	double value()const;

	template<class Archive>
	void save(Archive & ar, const unsigned int version) const;
	template<class Archive>
	void load(Archive & ar, const unsigned int version);
//
};

#endif /* SRC_INFNANFLOAT_H_ */
