/*!
 *  \brief      CombustionReactor modelling
 *  \details    This model assumes the CombustionReactor to be complete.
 *              The class has a template parameter `Delta` that can be:
 *               - the class `balanceMach` if only CombustionReactor is wanted
 *               - the class `Friction` if the friction needs to be included in the model.
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

#ifndef COMBUSTIONREACTOR_H_
#define COMBUSTIONREACTOR_H_

#include "core/propModule.h"
#include "cantera/zerodim.h"
#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/kinetics/GasKinetics.h"
#include "cantera/kinetics/importKinetics.h"
#include "cantera/base/stringUtils.h"
//#include "cantera/base/ctexceptions.h"
//#include "cantera/base/AnyMap.h"
//#include "cantera/IdealGasMix.h"

namespace hypro {
	class CombustionReactor : public propModule {
		std::tuple<double, double, vector<double>, double>
		plugFlow_reactor(double inlet_temp, double inlet_press, Cantera::compositionMap inlet_composition,
						 double inlet_massflow);

		std::tuple<double, double, Cantera::compositionMap>
		well_stirred_reactor(double inlet_temp, double inlet_press, Cantera::compositionMap inlet_composition,
							 double inlet_massflow, double reactor_volume);

	public:
		double length; ///<Length of the combustion reactor

		CombustionReactor();

		CombustionReactor(std::string name, Node *N1, Node *N2, int ID = defID_);

		CombustionReactor(const CombustionReactor &mod);

		virtual ~CombustionReactor();


		virtual double calculate();


		propModuleSERIAL(CombustionReactor);

		// virtual void serialize(std::ostream &os) const;

		// template<class Archive>
		// void serialize(Archive &ar, const unsigned int version);

	};
}
// BOOST_CLASS_EXPORT_KEY(hypro::CombustionReactor);

#endif /* COMBUSTIONREACTOR_H_ */
