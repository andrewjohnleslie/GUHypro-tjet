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

#ifndef ISENTROPICDUCT_H_
#define ISENTROPICDUCT_H_

#include "core/propModule.h"


//#define BOOST_MATH_INSTRUMENT
#include <boost/math/tools/roots.hpp>


namespace hypro {
class isentropicDuct: public propModule {
	std::unique_ptr<thermoKineticState> N1mod_; ///<N1 state modified considering pressure and temperature drop.

protected:
	double p02_p01_; ///<Stagnation pressure ratio between node `N1_` and node `N2_`
	double T02_T01_; ///<Stagnation temperature ratio between node `N1_` and node `N2_`

	/// Define the stagnation pressure loss across the duct
	/**Set the value of the pressure loss fraction `p02_p01_`
	 * This function needs to be reimplemented by subclasses.
	 * By default it sets the pressure loss fraction to 1.
	 */
	virtual void pressLoss();

	/// Define the stagnation temperature loss across the duct
	/**Set the value of the temperature loss fraction `T02_T01_`
	 * This function needs to be reimplemented by subclasses.
	 * By default it sets the temperature loss fraction to 1.
	 */
	virtual void thermLoss();

	///Modifies the N1mod_ state taking into account the pressure and thermal losses
	void applyLoss();

	///Returns the N1mod_ state
	const thermoKineticState& N1mod()const;

	///Returns the A2_A1 ratio of an isentropic expansion
	/**M2 is the Mach number at the end of isentropic expansion
	 * The isentropic law of the gas with variable Cp are implemented.
	 */
	double mfrEq(double& M2)const;

	///Returns the difference between A2_A1 ratio of an isentropic expansion and actual A2_A1 ratio of the module.
	/**M2 is the Mach number at the end of isentropic expansion
	 * No parameters are needed.
	 * This method can be used with class nonLinear to find M2.
	 */
	double mfrEq(double& M2, const void* par[]);

	///Function used to find the throat pressure.
	/**p is the tentative throat pressure
	 * No parameters are needed.
	 * The function returns the difference between the throat speed and the speed of sound.
	 */
	double pThroat(double& p, const void* par[]);

	///Function to find the pressure at node 2
	/**p is the tentative pressure
	 * No parameters are needed.
	 * The function returns the difference between the mfr of node 2 and the one of node 1.
	 */
	double mfrEqP(double& p2, const void* par[]);

public:
	///Set if the duct is supposed to have a sonic throat or not.
	/**If true the solution of the equation is searched assuming that
	 * `N2_->M() > 1` if `N2_->M() < 1` and `N2_->M() < 1` if `N2_->M() < 1`
	 * and vice-versa.
	 */
	bool choked_;
	bool solver_verbosity = false;

	isentropicDuct(){};  ///<Empty constructor
	isentropicDuct(std::string name, Node* N1, Node* N2, int ID=defID_);
	isentropicDuct(const isentropicDuct& mod); //copy constructor
	virtual ~isentropicDuct();

	virtual GPObject& duplicate () { return *(new isentropicDuct(*this)); }

	propModuleSERIAL(isentropicDuct);

	double calculate();

	///isentropic flow relations
	/**This relation is strictly valid for constant `gamma` gas only.
	 * It can be used for real gases accepting a certain error.
	 * In the actual version it is not used in this class, but it might be called from
	 * external classes.
	 *
	 * @param[in] Mach Mach number
	 * @param[in] gamma ratio of specific heats of the gas
	 * @param[out] p_p0 pressure static to stagnation ratio
	 * @param[out] T_T0 temperature static to stagnation ratio
	 */
	static void isentropic(
				const double& Mach, const double& gamma,
				double& p_p0, double& T_T0);

	///inverse of isentropic flow relations
	/**This relation is strictly valid for constant `gamma` gas only.
	 * It can be used for real gases accepting a certain error.
	 * In the actual version it is not used in this class, but it might be called from
	 * external classes.
	 *
	 * @param[in] gamma ratio of specific heats of the gas
	 * @param[in] p_p0 pressure static to stagnation ratio
	 * @return Mach number
	 */
	static double isentropicInv(const double& p_p0, const double& gamma);

	virtual void serialize(Archive& ar) const;
	virtual void unserialize(const Archive& ar);


};
}
//BOOST_CLASS_EXPORT_KEY(isentropicDuct);

#endif /* ISENTROPICDUCT_H_ */
