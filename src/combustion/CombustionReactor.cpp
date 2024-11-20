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

#include "CombustionReactor.h"
#include <cantera/base/Solution.h>
#include <vector>
#include <tuple>
#include <iostream>
#include <stdexcept>

namespace hypro {

    CombustionReactor::CombustionReactor(std::string name, Node *N1, Node *N2, int ID)
        : propModule(name, N1, N2, ID, 0) {}

    CombustionReactor::CombustionReactor(std::istream &is, const NodeMap &nodeMap, const ModuleMap &moduleMap)
        : propModule(is, nodeMap, moduleMap) {
        try {
            is >> length;
        } catch (std::exception &e) {
            throw std::runtime_error("Error: CombustionReactor input incorrect\nProvide Isp_ Me_ mfr_ mfrMax_.");
        }
    }

    CombustionReactor::CombustionReactor()
        : propModule(), length(0) {}

    CombustionReactor::CombustionReactor(const CombustionReactor &mod)
        : propModule(mod), length(mod.length) {}

    CombustionReactor::~CombustionReactor() {
        // Destructor
    }

    double CombustionReactor::calculate() {
        Node N1_inlet(*N1_);
        std::cout << N1_inlet.report() << std::endl;

        std::cout << "CombustionReactor\nN1 State:\n" << "Temperature: " << N1_inlet.getTemp() << ", Pressure: " << N1_inlet.getPress() << ", MFR: " << N1_inlet.mfr() << ", U: " << N1_inlet.getU() << std::endl;

        // Initialize Perfectly Stirred Reactor outputs
        double psr_pressure, psr_temperature;
        std::map<std::string, double> psr_composition;
        
        // Calculate volume of flame
        double reactor_volume = N1_inlet.A() * length * 0.01;
        std::tie(psr_pressure, psr_temperature, psr_composition) = well_stirred_reactor(N1_inlet.getTemp(),
                                                                                        N1_inlet.getPress(),
                                                                                        N1_inlet.getXMap(),
                                                                                        N1_inlet.mfr(), reactor_volume);

        std::cout << "Post-PSR State:\n" << "Temperature: " << psr_temperature << ", Pressure: " << psr_pressure << ", MFR: " << N1_inlet.mfr() << ", U: " << N1_inlet.getU() << std::endl;

        // Initialize Plugged-Flow Reactor outputs
        double pfr_pressure, pfr_temperature, pfr_speed;
        std::vector<double> pfr_composition;
        std::tie(pfr_pressure, pfr_temperature, pfr_composition, pfr_speed) = plugFlow_reactor(psr_temperature, psr_pressure,
                                                                                    psr_composition, N1_inlet.mfr() / N1_inlet.A());

        // Set Node 2 parameters
        N2_->setTPX(pfr_temperature, pfr_pressure, pfr_composition);
        N2_->setU(pfr_speed);
        std::cout << "N2 State:\n" << "Temperature: " << N2_->getTemp() << ", Pressure: " << N2_->getPress() << ", MFR: " << N2_->mfr() << ", U: " << N2_->getU() << std::endl;

        return 1.0;
    }

    std::tuple<double, double, std::vector<double>, double>
    CombustionReactor::plugFlow_reactor(double inlet_temp, double inlet_press,
                                        std::map<std::string, double> inlet_composition, double inlet_massflow) {

        bool verbosity = false;
        std::cout << "Started Plug Flow Reactor" << std::endl;
        
        //auto sol = Cantera::newSolution("gri30_highT.yaml", "gri30");
        auto sol = Cantera::newSolution("gri30_highT.yaml", "gri30", "none");
        auto gas = sol->thermo();
        auto kin = sol->kinetics();
        
        //sol->thermo()->setState_TPX(inlet_temp, inlet_press, inlet_composition);
	gas->setState_TPX(inlet_temp, inlet_press, inlet_composition);

        Cantera::FlowReactor pfr;
        double speed;

        // Dereference the shared_ptr to obtain a reference to the Kinetics object
        //pfr.setKineticsMgr(*sol->kinetics());
        pfr.setKineticsMgr(*kin);

        if (verbosity) {
            //std::cout << sol->thermo()->report() << std::endl;
            std::cout << gas->report() << std::endl;
        }

        pfr.setMassFlowRate(inlet_massflow / N1_->A());

        if (verbosity) {
            std::cout << pfr.energyEnabled() << std::endl;
        }

        Cantera::ReactorNet sim;
        sim.addReactor(pfr);

        double tnow = 0.0;

        if (verbosity) {
            std::cout << "Distance: " << pfr.distance() << std::endl;
            std::cout << "Speed: " << pfr.speed() << std::endl;
            std::cout << "Volume: " << pfr.volume() << std::endl;
        }

        while (pfr.distance() < length) {
            tnow = sim.step();
            std::cout << "tnow: " << tnow << ", distance: " << pfr.distance() << ", temperature: "
                      << pfr.temperature() << ", speed: " << pfr.speed() << std::endl;
            speed = pfr.speed();
        }

        if (verbosity) {
            //std::cout << sol->thermo()->report() << std::endl;
            std::cout << gas->report() << std::endl;
            std::cout << "Finished" << std::endl;
        }

        //std::vector<double> moles(sol->thermo()->nSpecies());
        std::vector<double> moles(gas->nSpecies());
        //sol->thermo()->getMoleFractions(moles.data());
        gas->getMoleFractions(moles.data());

        //return std::make_tuple(sol->thermo()->pressure(), sol->thermo()->temperature(), moles, speed);
        return std::make_tuple(gas->pressure(), gas->temperature(), moles, speed);
       
    }

    std::tuple<double, double, std::map<std::string, double>>
    CombustionReactor::well_stirred_reactor(double inlet_temp, double inlet_press,
                                            std::map<std::string, double> inlet_composition, double inlet_massflow,
                                            double reactor_volume) {
        bool verbosity = false;
        std::cout << "Started Perfectly Stirred Reactor" << std::endl;

        //auto sol = Cantera::newSolution("gri30_highT.yaml", "gri30");
        auto sol = Cantera::newSolution("gri30_highT.yaml", "gri30", "none");
        auto gas = sol->thermo();
        auto kin = sol->kinetics();
        
        //sol->thermo()->setState_TPX(inlet_temp, inlet_press, inlet_composition);
        gas->setState_TPX(inlet_temp, inlet_press, inlet_composition);

        Cantera::Reservoir fuel_in;
        //fuel_in.insert(*gas);	// Dereference the shared_ptr
        fuel_in.insert(sol);
        

        if (verbosity) {
            std::cout << sol->thermo()->report() << std::endl;
            //std::cout << gas->report() << std::endl;
        }

        Cantera::Reservoir exhaust;
        //exhaust.insert(*gas);	// Dereference the shared_ptr
        exhaust.insert(sol);

        Cantera::IdealGasReactor cstr;
        cstr.setInitialVolume(reactor_volume);
        //sol->thermo()->setState_TP(2000.00, inlet_press);
        //sol->thermo()->equilibrate("HP");
        gas->setState_TP(2000.00, inlet_press);
        gas->equilibrate("HP");

        if (verbosity) {
            //std::cout << sol->thermo()->report() << std::endl;
            std::cout << gas->report() << std::endl;
        }
        
        //cstr.insert(*gas);	// Dereference the shared_ptr
        cstr.insert(sol);

        // Dereference the shared_ptr to obtain a reference to the Kinetics object
        cstr.setKineticsMgr(*sol->kinetics());
        //cstr.setKineticsMgr(*kin());

        Cantera::MassFlowController m1;
        m1.install(fuel_in, cstr);
        m1.setMassFlowRate(inlet_massflow);

        Cantera::Valve v;
        v.install(cstr, exhaust);
        double Kv = 1.0;
        v.setValveCoeff(Kv);

        Cantera::ReactorNet sim;
        sim.addReactor(cstr);

        double tfinal = 1.0;
        double tnow = 0.0;
        double tres;

        while (tnow < tfinal) {
            tnow += 0.005;
            sim.advance(tnow);
            tres = cstr.mass() / v.massFlowRate();
            if (verbosity) {
                std::cout << tres << std::endl;
                std::cout << tnow << ", " << cstr.temperature() << std::endl;
            }
        }

        if (verbosity) {
            std::cout << exhaust.contents().report() << std::endl;
        }

        return std::make_tuple(exhaust.contents().pressure(), exhaust.contents().temperature(), 
                               exhaust.contents().getMoleFractionsByName());
    }
}


