/*!
 *  \brief      Isentropic variable area duct model
 *  \details    This class models a variable area duct under the assumption of isentropic flow.
 *              Non reversible effects can be modelled reimplementing the methods
 *              `pressLoss` and `thermLoss`.
 *              The flow is always considered at fixed composition.
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

#ifndef POLYTROPICDUCT_H_
#define POLYTROPICDUCT_H_

#include "core/propModule.h"

//#define BOOST_MATH_INSTRUMENT
#include <boost/math/tools/roots.hpp>


namespace hypro {
class polytropicDuct: public propModule {

	std::unique_ptr<thermoKineticState> N1mod_; ///<N1 state modified considering pressure and temperature drop.

protected:

	double p02_p01_; ///<Stagnation pressure ratio between node `N1_` and node `N2_`
	double T02_T01_; ///<Stagnation temperature ratio between node `N1_` and node `N2_`
	double eff_poly_; // polytropic efficiency for use in calculation


	///Returns the N1mod_ state
	const thermoKineticState& N1mod()const;


public:

	bool choked_ {false};
	bool solver_verbosity = false;

	polytropicDuct(){};  ///<Empty constructor
	polytropicDuct(std::string name, Node* N1, Node* N2, int ID=defID_);
	polytropicDuct(const polytropicDuct& mod); //copy constructor
	virtual ~polytropicDuct();

	virtual GPObject& duplicate () { return *(new polytropicDuct(*this)); }

	propModuleSERIAL(polytropicDuct);

	double calculate();


	double PmodeC(double& T2, const void *pass[]);
	double HmodeT(double& T2, const void *pass[]);

	double MfrEqST(double& M2, const void *pass[]);


	virtual void serialize(Archive& ar) const;
	virtual void unserialize(const Archive& ar);


};
}
//BOOST_CLASS_EXPORT_KEY(isentropicDuct);

#endif /* POLYTROPICDUCT_H_ */
