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

#ifndef INJECTIONPLATEPRESSURE_H_
#define INJECTIONPLATEPRESSURE_H_

#include "core/propModule.h"

namespace hypro {
/**Injection plate of a liquid propellant rocket.
 * It does not need an input node, the only inputs required are
 * the mass flow rate, the composition, the static temperature and the pressure injected mixture.
 * The velocity of the flow is calculated from the mfr and injected pressure.
 */
	class InjectionPlatePressure : public propModule {
		double mfr_; ///<Mass flow rate, can be reduced in case of feedBack
	public:
		double mfrMax_; ///<Design mass flow rate, i.e. the maximum possible
		vector<double> Y_; ///<Composition of injected mixture [Mass Fraction]
		double T_; ///<Static temperature of injected mixture
		double p_; ///<Pressure of the injected mixture

		InjectionPlatePressure();

		InjectionPlatePressure(std::string name, Node *N1, Node *N2, int ID = defID_);

		InjectionPlatePressure(const InjectionPlatePressure &mod);

		virtual ~InjectionPlatePressure();

		propModuleSERIAL(InjectionPlatePressure);

		virtual double calculate();

		vector<double> propellantMfr() const;

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
#endif /* INJECTIONPLATEPRESSURE_H_ */
