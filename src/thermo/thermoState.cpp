/*!
 *  \author     Robert Garner
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

#include "Archive.h"
#include "thermoState.h"
#include <boost/math/tools/roots.hpp>

#include "cantera/base/Solution.h"
#include "cantera/thermo/ThermoPhase.h"
#include <memory>
#include <string>
#include "numerics/nonLinear.h"
#include <iostream>

//#include <cassert>

//#include <boost/serialization/vector.hpp>

#define DO_QUOTE(X) #X
#define QUOTE(X) DO_QUOTE(X)

namespace hypro
{
bool thermoState::warning_ = false;
double thermoState::toll_ = 1e-6;


// Constructor to initialise the thermoState with a Cantera phase from a file
thermoState::thermoState(const std::string &infile, std::string id_) {
    auto sol = Cantera::newSolution(infile, id_);	// Creates a Solution object from the specified input file and phase id
    gasmodel_ = sol->thermo();			// Obtains the ThermoPhase (thermo()) object from the Solution
    infile_ = infile;
    equilibrium = false;
    
    std::cout << "thermoState created with infile: " << infile << " and id: " << id_ << std::endl;
}

// Method to create or reset the state with a new Cantera phase from a file
void thermoState::createState(const std::string &infile, std::string id_) {
    auto sol = Cantera::newSolution(infile, id_);	//
    gasmodel_ = sol->thermo();			//
    infile_ = infile;
    
    std::cout << "createState called with infile: " << infile << " and id: " << id_ << std::endl;
}

// Constructor to initialise the thermoState from an input stream
thermoState::thermoState(std::istream &is, const std::string &infile, std::string id_) {
    auto sol = Cantera::newSolution(infile, id_);	//
    gasmodel_ = sol->thermo();			//
    infile_ = infile;
    equilibrium = false;
    
    std::cout << "thermoState created from input stream with infile: " << infile << " and id: " << id_ << std::endl;
}

// Default constructor
thermoState::thermoState()
{
    std::cout << "No Cantera File Included";
}

// Copy constructor
thermoState::thermoState(const thermoState &ts) {
    auto sol = Cantera::newSolution(ts.infile_, ts.gasmodel_->name()); // Creates a Solution object from the input file and phase id of the source object
    gasmodel_ = sol->thermo();					  //
    infile_ = ts.getCTInfile();
    std::vector<double> state;					  // Saves & restores the state from the source object to this object
    ts.gasmodel_->saveState(state);
    gasmodel_->restoreState(state);

}

// Default destructor
thermoState::~thermoState() {
}

// Assignment operator
thermoState &thermoState::operator=(thermoState ts) {
    std::vector<double> state;			// Saves the state of the source object
    ts.gasmodel_->saveState(state);
    auto sol = Cantera::newSolution(ts.infile_, ts.gasmodel_->name());  //
    gasmodel_ = sol->thermo();					   //
    gasmodel_->restoreState(state);			// Restores the saved state to this object

    return *this;
}

// Gets the Cantera phase ID
std::string thermoState::getCTid() const {
    return gasmodel_->name();
}

// Gets the Cantera input file
std::string thermoState::getCTInfile() const {
    return infile_;
}

// Gets the temperature
double thermoState::getTemp() const {
    return gasmodel_->temperature();
}

// Gets the pressure
double thermoState::getPress() const {
    return gasmodel_->pressure();
}

// Sets the temperature and maintains equilibrium - IF NECESSARY!!!
void thermoState::setTemp(double newTemp) {
    if (!std::isfinite(newTemp) || newTemp <= 0) {
    	std::cout << "setTemp " << std::endl;
        throw std::invalid_argument("setTemp: Temperature must be a finite number greater than zero.");
    }

    gasmodel_->setState_TP(newTemp, getPress());
    if (equilibrium) {
        if (printflag) {			// If equilibrium is required, equilibrate the phase
            std::cout << "Equilibrate::setTemp " << std::endl;
        }
        try {
            gasmodel_->equilibrate("TP");
        } catch (Cantera::CanteraError& err) {
            fmt::print("TS::setTemp: {}\n", err.what());
        }
    }
    if (printflag) {
        std::cout << "setTemp: " << newTemp << " " << getTemp() << std::endl;
    }
    
//    update();
}

// Sets the pressure and maintains equilibrium - IF NECESSARY!!!
void thermoState::setPress(double newPressure) {
    if (!std::isfinite(newPressure) || newPressure <= 0) {
    	std::cout << "setPress " << std::endl;
        throw std::invalid_argument("setPress: Pressure must be a finite number greater than zero.");
    }

    // Sets the pressure while keeping the temperature constant
    gasmodel_->setState_TP(getTemp(), newPressure);
    if (equilibrium) {
        if (printflag) {
            std::cout << "Equilibrate::setPress " << std::endl;
        }
        try {
            gasmodel_->equilibrate("TP", "vcs");
        } catch (Cantera::CanteraError& err) {
            fmt::print("TS::setPress: {}\n", err.what());
        }
        
    }
    if (printflag) {
        std::cout << "setPress: " << newPressure << " " << getPress() << std::endl;
    }
   
//    update();
}

// Sets the temperature, pressure, and mole fractions
void thermoState::setTP(double newTemp, double newPressure) {
   if (!std::isfinite(newTemp) || newTemp <= 0) {
    	std::cout << "setTP - Temp " << std::endl;
        throw std::invalid_argument("setTP: Temperature must be a finite number greater than zero.");
    }
    if (!std::isfinite(newPressure) || newPressure <= 0) {
    	std::cout << "setTP - Press " << std::endl;
        throw std::invalid_argument("setTP: Pressure must be a finite number greater than zero.");
    }

    //std::cout << "Setting TP: Temperature = " << newTemp << ", Pressure = " << newPressure << std::endl;

    gasmodel_->setState_TP(newTemp, newPressure);
    if (equilibrium) {
        if (printflag) {
            std::cout << "Equilibrate::setTP " << std::endl;
        }
        try {
            gasmodel_->equilibrate("TP");
        } catch (Cantera::CanteraError& err) {
            fmt::print("TS::setTP: {}\n", err.what());
        }
    }
    if (printflag) {
        std::cout << "setTP: " << newTemp << " " << getTemp() << ", " << newPressure << " " << getPress() << std::endl;
    }
    
//    update();
}

// Sets the temperature, pressure, and mole fractions
void thermoState::setTPX(double newTemp, double newPressure, Cantera::compositionMap MoleFractionMap) {
    if (!std::isfinite(newTemp) || newTemp <= 0) {
    	std::cout << "setTPX1 - Temp " << std::endl;
        throw std::invalid_argument("setTPX: Temperature must be a finite number greater than zero.");
    }
    if (!std::isfinite(newPressure) || newPressure <= 0) {
    	std::cout << "setTPX1 - Pres " << std::endl;
        throw std::invalid_argument("setTPX: Pressure must be a finite number greater than zero.");
    }

    //std::cout << "Setting TPX: Temperature = " << newTemp << ", Pressure = " << newPressure << std::endl;
    
    gasmodel_->setState_TPX(newTemp, newPressure, MoleFractionMap);
//    update();
}

// Sets the temperature, pressure, and mole fractions
void thermoState::setTPX(double newTemp, double newPressure, std::vector<double> MoleFractionVector) {
    if (!std::isfinite(newTemp) || newTemp <= 0) {
    	std::cout << "setTPX2 - Temp " << std::endl;
        throw std::invalid_argument("setTPX: Temperature must be a finite number greater than zero.");
    }
    if (!std::isfinite(newPressure) || newPressure <= 0) {
    	std::cout << "setTPX2 - Pres " << std::endl;
        throw std::invalid_argument("setTPX: Pressure must be a finite number greater than zero.");
    }

    //std::cout << "Setting TPX: Temperature = " << newTemp << ", Pressure = " << newPressure << std::endl;
    
    gasmodel_->setState_TPX(newTemp, newPressure, MoleFractionVector.data());
//    update();
}

// Sets the temperature, pressure, and mole fractions
void thermoState::setTPY(double newTemp, double newPressure, std::vector<double> MassFractionVector) {
    if (!std::isfinite(newTemp) || newTemp <= 0) {
    	std::cout << "setTPY1 - Temp " << std::endl;
        throw std::invalid_argument("setTPY: Temperature must be a finite number greater than zero.");
    }
    if (!std::isfinite(newPressure) || newPressure <= 0) {
    	std::cout << "setTPY1 - Pres " << std::endl;
        throw std::invalid_argument("setTPY: Pressure must be a finite number greater than zero.");
    }

    //std::cout << "Setting TPY: Temperature = " << newTemp << ", Pressure = " << newPressure << std::endl;
    
    gasmodel_->setState_TPY(newTemp, newPressure, MassFractionVector.data());

//    update();
}

// Sets the temperature, pressure, and mole fractions
void thermoState::setTPY(double newTemp, double newPressure, Cantera::compositionMap MassFractionMap) {
    if (!std::isfinite(newTemp) || newTemp <= 0) {
    	std::cout << "setTPY2 - Temp " << std::endl;
        throw std::invalid_argument("setTPY: Temperature must be a finite number greater than zero.");
    }
    if (!std::isfinite(newPressure) || newPressure <= 0) {
    	std::cout << "setTPY2 - Pres " << std::endl;
        throw std::invalid_argument("setTPY: Pressure must be a finite number greater than zero.");
    }

    //std::cout << "Setting TPY: Temperature = " << newTemp << ", Pressure = " << newPressure << std::endl;
    
    gasmodel_->setState_TPY(newTemp, newPressure, MassFractionMap);

//    update();
}


// Sets the enthalpy and pressure and maintains equilibrium - IF NECESSARY!!!
void thermoState::setHP(double newH, double newPressure) {
    
    gasmodel_->setState_HP(newH, newPressure);
    if (equilibrium) {
        if (printflag) {
            std::cout << "Equilibrate::setHP " << std::endl;
        }
        try {
            gasmodel_->equilibrate("HP");
        } catch (Cantera::CanteraError& err) {
            fmt::print("TS::setHP: {}\n", err.what());
        }
    }

    //std::cout << "Setting HP: Enthalpy = " << newH << ", Pressure = " << newPressure << std::endl;    

//    update();
}
 
// Sets the enthalpy and entropy and maintains equilibrium - IF NECESSARY!!!
void thermoState::setSH(double S, double H) {
    if (!std::isfinite(S) || !std::isfinite(H)) {
    	std::cout << "setSH - Entr " << std::endl;
        throw std::invalid_argument("setSH: Entropy or Enthalpy is not a finite number.");
    }
    
    //std::cout << "Setting SH: Entropy = " << S << ", Enthalpy = " << H << std::endl;

    double guess = 100000;
    double factor = 2;

    const boost::uintmax_t maxit = 30;
    boost::uintmax_t it = maxit;
    bool is_rising = true;
    int digits = std::numeric_limits<double>::digits;
    int get_digits = digits - 3;
    boost::math::tools::eps_tolerance<double> tol(get_digits);

    // Residual function: error in enthalpy as a function of Pressure
    auto h_err = [this, S, H](double P) {
        gasmodel_->setState_SP(S, P);
        return h() - H;
    };

    std::pair<double, double> r = boost::math::tools::bracket_and_solve_root(h_err, guess, factor, is_rising, tol,
                                                                             it);
    double P_calc = r.first + (r.second - r.first) / 2;

    gasmodel_->setState_SP(S, P_calc);
    if (equilibrium)
    {
        if (printflag)
        {
            std::cout << "Equilibrate::setSH " << std::endl;
        }
        gasmodel_->equilibrate("SP");
    }

    if (printflag)
    {
        std::cout << "setSH: " << S << " " << s() << "," << H << " " << h() << std::endl;
    }
}

// Sets the entropy and temperature
void thermoState::setST(double S, double T) {
    if (!std::isfinite(S)) {
    	std::cout << "setST - Entr " << std::endl;
        throw std::invalid_argument("setST: Entropy is not a finite number.");
    }
    /* if (!std::isfinite(T) || T <= 0) {
    	std::cout << "setST - Temp " << std::endl;
        throw std::invalid_argument("setST: Temperature must be a finite number greater than zero.");
    }
    */

    //std::cout << "Setting ST: Entropy = " << S << ", Temperature = " << T << std::endl;

    double guess = 100;
    double factor = 2;

    const boost::uintmax_t maxit = 30;
    boost::uintmax_t it = maxit;
    bool is_rising = true;
    int digits = std::numeric_limits<double>::digits;
    int get_digits = digits - 3;
    boost::math::tools::eps_tolerance<double> tol(get_digits);

    // Residual function: error in Temperature as a function of Pressure
    auto T_err = [this, S, T](double P) {
        gasmodel_->setState_SP(S, P);
        if (equilibrium) {
            try {
                gasmodel_->equilibrate("SP");
            } catch (Cantera::CanteraError& err) {
                fmt::print("TS::setST: {}\n", err.what());
            }
        }
        return getTemp() - T;
    };

    // Finds the pressure that balances the temperature using the root solver
    std::pair<double, double> r = boost::math::tools::bracket_and_solve_root(T_err, guess, factor, is_rising, tol, it);
    double newP = r.first + (r.second - r.first) / 2;

    if (printflag) {
        std::cout << "setST: " << S << " " << s() << "," << T << " " << getTemp() << std::endl;
    }

    setTP(T, newP);		// Sets the state using the calculated pressure
}

// Sets the entropy and pressure
void thermoState::setSP(double S, double P) {
    if (!std::isfinite(S)) {
    std::cout << "setSP - Entr " << std::endl;
        throw std::invalid_argument("setSP: Entropy is not a finite number.");
    }
    if (!std::isfinite(P) || P <= 0)  {
    	std::cout << "setSP - Pres " << std::endl;
        throw std::invalid_argument("setSP: Pressure must be greater than zero.");
    }

    //std::cout << "Setting SP: Entropy = " << S << ", Pressure = " << P << std::endl;

    gasmodel_->setState_SP(S, P);
    if (equilibrium) {
        try {
            gasmodel_->equilibrate("SP");
        } catch (Cantera::CanteraError& err) {
            fmt::print("TS::setSP: {}\n", err.what());
        }
    }
    if (printflag) {
        std::cout << "setSP: " << S << " " << s() << "," << P << " " << getPress() << std::endl;
    }
//    update();
}

// Sets the phase to equilibrium
void thermoState::set_to_equilibrium() {
    gasmodel_->equilibrate("HP");
}

// Gets the density
double thermoState::rho() const {
    return gasmodel_->density();
}

// Gets the enthalpy
double thermoState::h() const {
    if (gasmodel_->temperature() > 6000)
    {
        if (warning_)
        {
            std::cout << "Temperature out of range of data: " << gasmodel_->temperature() << " K" << std::endl;
        }
    }
    return gasmodel_->enthalpy_mass();
}


double thermoState::Cp() const {
    if (gasmodel_->temperature() > 6000)
    {
        if (warning_)
        {
            std::cout << "Temperature out of range of data: " << gasmodel_->temperature() << " K" << std::endl;
        }
    }
    return gasmodel_->cp_mass();
}

// Gets the ratio of specific heats (gamma)
double thermoState::gamma() const {
    return gasmodel_->cp_mole() / gasmodel_->cv_mole();
}

// Gets the speed of sound
double thermoState::geta() const {
    return sqrt(gamma() * gasmodel_->RT() / gasmodel_->meanMolecularWeight());
}

// Gets the product of gas constant and temperature
double thermoState::RT() const {
    return gasmodel_->RT();
}

// Gets the molecular weight
double thermoState::W() const {
    return gasmodel_->meanMolecularWeight();
}

// Gets the gas constant (R)
double thermoState::R() const {
    return Cantera::GasConstant;
}

// Gets the entropy per unit mass
double thermoState::s() const {
    return gasmodel_->entropy_mass();
}

// Gets the enthalpy
double thermoState::getEnthalpy() const {
    return gasmodel_->enthalpy_mass();
}

/**
 * 	Printing the molecularweights of all the species available in the simulation
 *	for (uint i=0; i<gas.WX().size(); i++)
 *	cout << gas.WX()[i] << endl;
 *	cout << gas.WX().size();
 *
 */
// Gets the molecular weights of all species
std::vector<double> thermoState::WX() const {
    return gasmodel_->molecularWeights();
}

/**
 *	for (uint i=0; i<gas.X().size(); i++)
 *	cout << gas.X()[i] << " " << gas.WX()[i] << endl;
 *	cout << gas.X().size() << endl;
 */
// Gets the names of all species
std::vector<std::string> thermoState::XbyNames() const {
    return gasmodel_->speciesNames();
}

// Gets the mole fractions of all species
std::vector<double> thermoState::X() const {
    std::vector<double> moles(gasmodel_->nSpecies());
    gasmodel_->getMoleFractions(&moles[0]);
    return moles;
}

// Sets the mole fractions of all species
void thermoState::X(std::vector<double> mole_fractions) {
    gasmodel_->setState_TPX(gasmodel_->temperature(), gasmodel_->pressure(), mole_fractions.data());
}

// Gets the mole fraction of a specific species
double thermoState::X(std::string speciesName) {
    Cantera::compositionMap speciesMap = gasmodel_->getMoleFractionsByName(0);
    return speciesMap[speciesName];
}

// Gets the mass fraction of a specific species
double thermoState::Y(std::string speciesName) {
   return gasmodel_->getMassFractionsByName(0)[speciesName];
}

// Gets the mole fractions as a composition map
Cantera::compositionMap thermoState::getXMap() const {
    return gasmodel_->getMoleFractionsByName(0);
}

// Gets the mass fractions as a composition map
Cantera::compositionMap thermoState::getYMap() const {
    return gasmodel_->getMassFractionsByName(0);
}

/**
	 * Function calculates the temperature under isentropic conditions
	 * from one state to another.
	 *
	 * Example of use:
	 * ThermoKineticEquiv gas("gri30.cti", "gri30");
	 * gas.setState_TP(1000,101325);
	 * ThermoKineticEquiv gas2("gri30.cti", "gri30");
	 * gas2.setState_TP(2000,201325);
	 * gas2.isentropicP(gas, 201325);
	 * Results in gas2 with a pressure of 201325, an entropy equivalent to that of gas1
	 * and an altered temperature.
	 */
// Sets the state under isentropic conditions knowing the final pressure	 
void thermoState::isentropicP(const thermoState &ST1, double p) {
    try
    {
        gasmodel_->setState_SP(ST1.s(), p);
    }
    catch (const Cantera::CanteraError& e) // Catching by reference
    {
        std::cout << "Error in isentropicP: " << e.what() << std::endl; // Include error message
    }
  
//    update();
}

/**
	 * Function calculates the pressure under isentropic conditions from
	 * one state to another state. Input parameters is the object representing state 1
	 * and the desired temperature.
	 *
	 * Example of use:
	 * ThermoKineticEquiv gas("gri30.cti", "gri30");
	 * gas.setState_TP(1000,101325);
	 * ThermoKineticEquiv gas2("gri30.cti", "gri30");
	 * gas2.setState_TP(2000,201325);
	 * gas2.isentropicT(gas, 2000);
	 * Results in gas2 with a temperature of 2000, an entropy equivalent to that of gas1
	 * and an altered pressure.
	 */
// Sets the state under isentropic conditions knowing the final temperature	 
void thermoState::isentropicT(const thermoState &ST1, double T) {
    double S = ST1.s();
    const void *par[2] = {&S, &T};
    
    // Uses a nonlinear solver to find the pressure that maintains entropy
    double p = nonLinear<thermoState>::funSolve(&thermoState::s, 1, 1e+8, *this, par);
    gasmodel_->setState_TP(T, p);
    if (equilibrium) {
        gasmodel_->equilibrate("TP");
    }
//    update();
}

/**
 *	Method to calculate the entropy used in the nonlinear function solver
 *	when the isentropicT is called.
 */
// Calculates the entropy for use in nonlinear function solver 
double thermoState::s(double &P, const void *par[]) {
    double &S = *(double *)par[0];
    double &T = *(double *)par[1];
    gasmodel_->setState_TP(T, P);
    if (equilibrium){
        gasmodel_->equilibrate("TP");
    }
    return s() - S;
}

// Checks if the temperature is within the valid range
void thermoState::limit(double T) const {
    if (getTemp() < gasmodel_->minTemp() || getTemp() > gasmodel_->maxTemp()) {
        if (warning_)
        {
            std::cout << "Attempt to use thermodynamic model out of range:" << gasmodel_->maxTemp() << " -> " << gasmodel_->maxTemp() << "@" << getTemp() << std::endl;
        }
    }
}

// Checks if the temperature and pressure are positive
inline void thermoState::check() const {
    if (getTemp() <= 0)
    {
        throw std::range_error("Error: Pressure lower than 0");
    }
    if (getPress() <= 0)
    {
        throw std::range_error("Error: Temperature lower than 0");
    }
}

// Gets a report of the state
std::string thermoState::report() const {
    return gasmodel_->report(true, -1e-14);
}

// Sets the state using enthalpy and pressure
void thermoState::Th(const double &h) {
    gasmodel_->setState_HP(h, getPress());
    //logGasModelState(gasmodel_);
//    update();
}

void thermoState::update() {
}

// Gets the lowest valid temperature
double thermoState::Tlow() const {
    return gasmodel_->minTemp();
}

// Gets the highest valid temperature
double thermoState::Thigh() const {
    return gasmodel_->maxTemp();
}

// Serialise state
void thermoState::serialize(Archive &ar) const {
    ar.put("Infile", getCTInfile());
    ar.put("PhaseID", getCTid());
    ar.put("Pressure", getPress());
    ar.put("Temperature", getTemp());
    ar.put("Composition", X());
}

// Unserialise state
void thermoState::unserialize(const Archive &ar)
{
    createState(ar.get<std::string>("Infile"), ar.get<std::string>("PhaseID"));
    gasmodel_->setState_TPX(ar.getDouble("Temperature"), ar.getDouble("Pressure"), ar.getVector("Composition").data());
    //logGasModelState(gasmodel_);
}
}



