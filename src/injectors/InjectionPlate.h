/*!
 *  \brief      Injection plate of a liquid propellant rocket.
 *  \details     It does not need an input node, the only inputs required are
 *              the mass flow rate, the composition, the static temperature and the Mach number of the injected mixture.
 *              The total and static pressure are calculated based on the mfr.
 *              The Mach number default value is 0.1 when the injection is reduced due to a choking feedback
 *              the Mach is kept constant, indeed it depends solely by the ratio between the chamber area and the
 *              throat area that is usually constant.
 *              In other words changing the Mach number is equivalent to changing this ratio.
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

#ifndef INJECTIONPLATE_H_
#define INJECTIONPLATE_H_

#include "core/propModule.h"

namespace hypro {
	class InjectionPlate : public propModule {
		double mfr_; ///<Mass flow rate, can be reduced in case of feedBack [kg/s]
	public:
		double mfrMax_; ///<Design mass flow rate, i.e. the maximum possible [kg/s]
		std::vector<double> Y_; ///<Composition of injected mixture [Mass Fraction]
		double T_; ///<Static temperature of injected mixture [K]
		double M_; ///<Mach number of the injected mixture

		InjectionPlate();

		InjectionPlate(std::string name, Node *N1, Node *N2, int ID = defID_);

		InjectionPlate(const InjectionPlate &mod);

		virtual ~InjectionPlate();

		propModuleSERIAL(InjectionPlate);

		virtual GPObject &duplicate() { return *(new InjectionPlate(*this)); }

		virtual double calculate();

		std::vector<double> propellantMfr() const;

		void reduce(const double &mfr);

		double reduce() const;

		bool isreduceable() const;

		void unreduce();

	virtual void serialize(Archive& ar) const;
	virtual void unserialize(const Archive& ar);

	virtual bool canChoke()const;

	virtual std::vector<const Node*> inNode()const;


};
}
//BOOST_CLASS_EXPORT_KEY(InjectionPlate);

#endif /* INJECTIONPLATE_H_ */
