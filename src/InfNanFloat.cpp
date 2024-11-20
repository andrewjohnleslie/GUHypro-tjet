/*
 * InfNanFloat.cpp
 *
 *  Created on: 8 Mar 2017
 *      Author: trb12187
 */

#include <InfNanFloat.h>
#include <math.h>
//#include <boost/serialization/nvp.hpp>
#include <limits>
#include <cmath>
//#include <boost/archive/xml_oarchive.hpp>
//#include <boost/archive/xml_iarchive.hpp>

InfNanFloat::InfNanFloat(double value):
	value_(value),
	isInf_(value_ == std::numeric_limits<double>::infinity()),
	isNan_(std::isnan(value_)){
}

InfNanFloat::~InfNanFloat() {
	// TODO Auto-generated destructor stub
}





double InfNanFloat::value()const{
	return value_;
}



