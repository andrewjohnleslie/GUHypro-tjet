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

#ifndef ROCKET_H_
#define ROCKET_H_

#include "core/propModule.h"

namespace hypro {
	class Rocket : public propModule {

		double mfr_; ///<Mass flow rate, can be reduced in case of feedBack
	public:
		double mfrMax_; ///<Design mass flow rate, i.e. the maximum possible
		double Isp_; ///<Specific impulse
		double Me_; ///<Nozzle exit Mach number

		Rocket();

		Rocket(std::string name, Node *N1, Node *N2, int ID = defID_);

		Rocket(const Rocket &mod);

		virtual ~Rocket();

		virtual GPObject &duplicate() { return *(new Rocket(*this)); }

		propModuleSERIAL(Rocket);

		virtual double calculate();

		std::vector<double> propellantMfr() const;

		void reduce(const double &mfr);

		double reduce() const;

		///Test whether the module is reducible
		/**@todo maybe fix the function spelling isreducible instead of isreduceable
         * @todo probably not applicable to this version of the code. Check this out.
         *
         * @return true if the module has all the reduce methods defined
         */
		bool isreduceable() const;

		void unreduce();

	virtual void serialize(Archive& ar) const;
	virtual void unserialize(const Archive& ar);

	virtual bool canChoke()const;

	virtual std::vector<const Node*> inNode()const;


};
}
//BOOST_CLASS_EXPORT_KEY(Rocket);

#endif /* ROCKET_H_ */
