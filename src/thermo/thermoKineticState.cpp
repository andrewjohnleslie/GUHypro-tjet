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
#include "thermoKineticState.h"
//#define BOOST_MATH_INSTRUMENT
#include <boost/math/tools/roots.hpp>


namespace hypro {
    thermoKineticState::thermoKineticState(const std::string &infile, std::string id_) :
            thermoState(infile, id_),
            U_(0),
            stagnation_(*this) {
    }

    thermoKineticState::thermoKineticState(std::istream &is, const std::string &infile, std::string id_) :
            thermoState(is, infile, id_),
            stagnation_(*this) {
        is >> U_;
    }

    thermoKineticState::thermoKineticState() {
        std::cout << "No Cantera File Included";
    }

    thermoKineticState::thermoKineticState(const thermoKineticState &tks) :
            thermoState(tks),
            U_(tks.U_),
            stagnation_(tks.stagnation_) {
    }

    thermoKineticState::~thermoKineticState() {

    }

    thermoKineticState &thermoKineticState::operator=(thermoKineticState ts) {
        thermoState::operator=(ts);
        U_ = ts.U_;

        return *this;
    }

    double thermoKineticState::M() const {
        return U_ / geta();
    }

    double thermoKineticState::G() const {
        return U_ * rho();
    }

    double thermoKineticState::I() const {
        return U_ * G() + getPress();
    }

    double thermoKineticState::H() const {
        return 0.5 * std::pow(U_, 2.0) + h();
    }

    double thermoKineticState::T0() const {
        return stagnation_.getTemp();
    }

    double thermoKineticState::getT0() {
        update();
        return stagnation_.getTemp();
    }

    double thermoKineticState::p0() const {
//        fmt::print(stagnation_.report());
        return stagnation_.getPress();
    }

    double thermoKineticState::getp0() {
        update();
        return stagnation_.getPress();
    }

    void thermoKineticState::setU(double newU_) {
        U_ = newU_;
//        update();
    }

    double thermoKineticState::getU() const {
        return U_;
    }

    void thermoKineticState::setEquilibrium(bool flag) {
        equilibrium = flag;
        //TODO added stagnation_.equilibrium back in
//        stagnation_.equilibrium = flag;
    }


    void thermoKineticState::update() {
        if (calculate_stagnation) {

            thermoState::update();
            stagnation_.setTPX(getTemp(), getPress(), X());

            if (printflag) {
                std::cout << "U= " << getU() << std::endl;
                std::cout << "Normal State: " << std::endl;
                std::cout << report() << std::endl;
            }

            updateStagnation();

            if (printflag) {
                std::cout << "Stagnation State: " << std::endl;
                std::cout << stagnation_.report() << std::endl;
                std::cout << "############################################################" << std::endl;
            }
        }

    }


    void thermoKineticState::updateStagnation() {
        double totalH = H();
        double S = s();

        if (printflag) {
            std::cout << "updateStagnation() -----------------------------------" << std::endl;
            std::cout << "totalH: " << totalH << " S: " << S << std::endl;
        }
        stagnation_.setSH(S, totalH);
    }

    void thermoKineticState::setStagnation(const thermoState &ST1, double M) {
        double H1 = ST1.h();
        double S1 = ST1.s();

        if (printflag) {
            std::cout << "setStagnation -----------------------------------" << std::endl;
            std::cout << "Stagnation Param: " << ST1.getTemp() << " " << ST1.getPress() << std::endl;
        }
        double guess = 101325;
        double factor = 2;

        const boost::uintmax_t maxit = 20;
        boost::uintmax_t it = maxit;
        bool is_rising = false;
        int digits = std::numeric_limits<double>::digits;
        int get_digits = digits - 3;
        boost::math::tools::eps_tolerance<double> tol(get_digits);

        // Residual function: error in Enthalpy as geta function of Temperature
        auto h_err = [this, M, H1, S1](double p) {
            try {
//                gasmodel_->setState_TP(T, getPress());
//                if (equilibrium) {
//                    gasmodel_->equilibrate("TP");
//                }
//                setST(S1, T);
                setSP(S1, p);

//                setTP(T, getPress());
            } catch (Cantera::CanteraError& err) {
                fmt::print("TKS: {}\n", err.what());
            }
            if (printflag) {
                fmt::print("M: {}, H1: {}, S1: {}, p: {}, T: {}, res: {}\n", M, H1, S1, p, getTemp(),
                           H1 - h() - (pow(M, 2) * pow(geta(), 2)) / 2);
            }
            return H1 - h() - (pow(M, 2) * pow(geta(), 2)) / 2;
        };

        std::pair<double, double> r = boost::math::tools::bracket_and_solve_root(h_err, guess, factor, is_rising, tol,
                                                                                 it);

        double newp = r.first + (r.second - r.first) / 2;

//        setST(S1, newT);
//        setST(S1, newT);
        setSP(S1, newp);
        if (printflag) {
            fmt::print("Selected p: {}\n", newp);
            fmt::print("M: {}, H1: {}, S1: {}, p: {}, T: {}, res: {}\n", M, H1, S1, getPress(), getTemp(),
                       H1 - h() - (pow(M, 2) * pow(geta(), 2)) / 2);
        }
        U_ = M * geta();
        if (printflag) {
            std::cout << "setStagnation: " << "H1: " << H1 << " h: " << h() << " S1: " << S1 << " s: " << s()
                      << " newT: " << newp << std::endl;
        }
    }

thermoState& thermoKineticState::getstag() {
    update();
    return stagnation_;
}

void thermoKineticState::serialize(Archive& ar) const{
	thermoState::serialize(ar);
	ar.put("StagnationTemperature", T0());
	ar.put("StagnationPressure", p0());
	ar.put("Enthalpy", H());
	ar.put("Entropy", s());
	ar.put("Mach" ,M());
	ar.put("Velocity", U_);

}

void thermoKineticState::unserialize(const Archive& ar){
	thermoState::unserialize(ar);

	U_ = ar.getDouble("Velocity");
}
}
