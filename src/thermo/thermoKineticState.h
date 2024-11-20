/*!
 *  \brief      Thermo-kinetic state
 *  \author     Alessandro Mogavero, Robert Garner
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

#ifndef HYPRO_THERMOKINETICSTATEC_H
#define HYPRO_THERMOKINETICSTATEC_H

#include "thermoState.h"
#include "numerics/nonLinear.h"
#include <iostream>


namespace hypro {
    class thermoKineticState : public thermoState {
        double U_; ///<Speed
        bool printflag = false;
    protected:
        thermoState stagnation_;
    public:
        thermoKineticState(const std::string &infile, std::string id_);

        thermoKineticState(std::istream &is, const std::string &infile, std::string id_);

        thermoKineticState();

        thermoKineticState(const thermoKineticState &tks);

        virtual ~thermoKineticState();

        //assignment operator
        thermoKineticState &operator=(thermoKineticState ts);

        ///Mach [-]
        double M() const;

        ///Mass flow [kg/(s*m^2)]
        double G() const;

        ///Impulse (i.e. momentum flow) [Pa]
        double I() const;

        ///Total enthalpy [J/kg]
        double H() const;

        ///Total temperature [K]
        double T0() const;

        double getT0();

        ///Total pressure [Pa]
        /**The value is calculated considering the entropy of an ideal gas.
         * Note however that no hypothesis is made on the specific heat relation.
         */
        double p0() const;

        double getp0();

        void setU(double newU_);

        double getU() const;

        virtual void update();

        void setEquilibrium(bool flag);

        void updateStagnation();
        //void serialize(std::ostream& os)const;

        void setStagnation(const thermoState &ST1, double M);

        //double findT(double& T, const void* par[]);
        thermoState &getstag();

        virtual void serialize(Archive& ar) const;
	    virtual void unserialize(const Archive&);
    };
}

#endif //HYPRO_THERMOKINETICSTATEC_H
