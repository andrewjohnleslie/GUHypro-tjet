/*!
 *  \brief      Thermal state
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

#ifndef THERMOSTATE_H_
#define THERMOSTATE_H_

#include "cantera/base/Solution.h"
//#include "cantera/thermo/ThermoPhase.h"
//#include "cantera/thermo/ThermoFactory.h"
//#include "cantera/thermo/Phase.h"
//#include "cantera/IdealGasMix.h"
#include <memory>
#include <string>
#include "numerics/nonLinear.h"
#include <iostream>

//#include <cassert>

namespace hypro
{
class Archive;
class thermoState
{


    /**It is an extension of the method S of the gas model,
         * that considers also variation of pressure and that works also out of T limits.
         *
         * @todo This would be better placed in the gas model, but major changes are needed.
         *
         * @param T temperature
         * @param p pressure
         */
    //Returns the entropy at given Temperature and pressure     
    double s(double &T, const void *par[]);
    
    // New helper function for logging gas model state
    //void logGasModelState(const std::shared_ptr<Cantera::ThermoPhase>& gasmodel_) const;

    bool printflag = false;
private:
	double lastH_ = 0.0;
protected:
   static double toll_; ///< Tolerance for numerical calculations.
    void limit(double T) const;
    void check() const;

public:
    static bool warning_; ///< Whether to display the warnings
    std::shared_ptr<Cantera::ThermoPhase> gasmodel_;  
    bool calculate_stagnation = true;
    bool equilibrium = false;
    std::string infile_;
    //double lastH_;

    thermoState();

    /// Constructor
    /**
     * @param specieGas gas model
     */
    thermoState(const std::string &infile, std::string id_);

    /// Construct from input stream
    /**
     * This constructor is protected since it is only meant to be called by
     * the unserialize method.
     *
     * @param is input stream
     * @param specieGas gas model
     */
    thermoState(std::istream &is, const std::string &infile, std::string id_);

    thermoState(const thermoState &ts);

    // Accessors
    std::string getCTid() const;
    std::string getCTInfile() const;
    
    // Destructor
    virtual ~thermoState();

    // Assignment operator
    thermoState &operator=(thermoState ts);

    // State creation
    void createState(const std::string &infile, std::string id_);

    // Get X_
    std::vector<std::string> XbyNames() const;	// Accessor
    /// Get X_
    std::vector<double> X() const;			// Accessor
    /// Set X_
    void X(std::vector<double>);
    double X(std::string);				// Accessor
    double Y(std::string);				// Accessor

    // Accessors
    double getTemp() const;
    double getPress() const;

    // Mutators
    void setTemp(double newTemp);
    void setPress(double newPressure);
    void setTP(double newTemp, double newPressure);
    void setTPX(double newTemp, double newPressure, Cantera::compositionMap MoleFractionMap);
    void setTPX(double newTemp, double newPressure, std::vector<double> MoleFractionVector);
    void setTPY(double newTemp, double newPressure, Cantera::compositionMap MoleFractionMap);
    void setTPY(double newTemp, double newPressure, std::vector<double> MoleFractionVector);
    //void setHP(double newH, double newPressure);
    void setSP(double S, double P);

    // Accessors
    Cantera::compositionMap getXMap() const;
    Cantera::compositionMap getYMap() const;

    // Mutators
    void setST(double S, double T);
    void setSH(double S, double H);
    void setHP(double newH, double newPressure);
    void set_to_equilibrium();

    // Accessors
    /// Return density [kg/m^3]
    double rho() const;
    /// Enthalpy [J/kg]
    double h() const;
    /// Heat capacity at constant pressure [J/(kg K)]
    double Cp() const;
    /// gamma = cp/cv []
    double gamma() const;
    /// Speed of sound [m/s]
    double geta() const;
    /// Gas constant R
    double R() const;
    /// Average molecular Weight
    double W() const;
    /// Entropy per unit mass
    double s() const;
    double RT() const;
    double getEnthalpy() const;
    
    // Accessors
    std::vector<double> WX() const;
    std::string report() const;

    /// Set state as coming from an isentropic transformation knowing the temperature
    /**
     * @param ST1 is the initial state.
     * @param T is the temperature at the end of transformation.
     */
    void isentropicT(const thermoState &ST1, double T);	// Transformation

    /// Set state as coming from an isentropic transformation knowing the pressure
    /**
     * @param ST1 is the initial state.
     * @param p is the pressure at the end of transformation.
     */
    void isentropicP(const thermoState &ST1, double p);	// Transformation

    /// Temperature from Enthalpy
    /** Calculates the temperature numerically
     *
     * @param h static enthalpy
     * @param T0 initial tentative temperature
     */
    void Th(const double &h);					// Transformation

    // Update
    virtual void update();

    // Serialisation
    void serialize(std::ostream &os) const;
    
    // Accessors
    double Tlow() const;
    double Thigh() const;

    /// Write to stream - Serialisation
    virtual void serialize(Archive &ar) const;
    virtual void unserialize(const Archive &ar);

    // Serialisation
    template <class Archive>
    void save(Archive &ar, const unsigned int version) const;
    template <class Archive>
    void load(Archive &ar, const unsigned int version);
    
};
}

#endif /* THERMOSTATE_H_ */
