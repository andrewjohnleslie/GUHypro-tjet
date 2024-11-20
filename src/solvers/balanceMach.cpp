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

#include "balanceMach.h"
#include "nonLinear.h"
#include <boost/math/tools/roots.hpp>


namespace hypro {
	balanceMach::balanceMach(std::string name, Node *N1, Node *N2, int ID)
			:
			propModule(name, N1, N2, ID),
			G1_(0),
			I1_(0),
			H1_(0) {
	}

	balanceMach::balanceMach(std::string name, Node *N1, Node *N2, int ID, int args)
			:
			propModule(name, N1, N2, ID, args),
			G1_(0),
			I1_(0),
			H1_(0) {
	}

	balanceMach::balanceMach()
			:
			propModule(),
			G1_(0),
			I1_(0),
			H1_(0) {
	}

	balanceMach::balanceMach(std::istream &is, const NodeMap &nodeMap, const ModuleMap &moduleMap)
			:
			propModule(is, nodeMap, moduleMap) {
	}

	balanceMach::balanceMach(const balanceMach &mod)
			:
			propModule(mod),
			G1_(mod.G1_),
			I1_(mod.I1_),
			H1_(mod.H1_) {
	}

	balanceMach::~balanceMach() {
		// TODO Auto-generated destructor stub
	}

	void balanceMach::checkAreas() {
		N2_->A(N1_->A());
		N2_->Amax_ = N2_->A();
		N2_->Amin_ = N2_->A();
	}

	double balanceMach::calculate() {
		checkAreas();

    // std::cout << chokingMach_ << std::endl;
		G1_ = N1_->G() * N1_->A() / N2_->A();
		I1_ = N1_->I() * N1_->A() / N2_->A();
		H1_ = G1_ * N1_->H() * N1_->A() / N2_->A();

		//Initialize node N2_
		N2_->setTP(N1_->getTemp(), N1_->getPress());
		N2_->setU(N1_->getU());
//        std::cout << N1_->Name_ << std::endl;
//        std::cout << N2_->Name_ << std::endl;
//        std::cout << "BalanceMach::calculate::changeComposition():" << N2_->report() << std::endl;
//		fmt::print(N1_->report());
//		fmt::print(N2_->report());
		changeComposition();
//		fmt::print(N1_->report());
//		fmt::print(N2_->report());

//        std::cout << "BalanceMach::calculate::changeComposition():" << N2_->report() << std::endl;

		double rChok = choked();
		if (solver_verbosity){
			fmt::print("\t{}\n", rChok);
		}

		if (rChok < 0) {
			return rChok;
		}

		double Mmin, Mmax, guess;
		if (N1_->M() <= 1) {
			Mmin = 0.0;
			Mmax = 1.0;
            guess = 0.1;
		} else {
//			Mmin = 1.0;
            Mmin = 1;
			Mmax = 20.0;
            guess = 1;
		}


        bool use_boost = true;
        double M2;

        if (use_boost == true){
           // std::cout << "Using Boost" << std::endl;
			if (solver_verbosity){
				fmt::print("Main balanceMach loop\n");
			}

            double factor = 1.25;

            const boost::uintmax_t maxit = 30;
            boost::uintmax_t it = maxit;
            bool is_rising = true;
            int digits = std::numeric_limits<double>::digits;
            int get_digits = digits - 3;
            boost::math::tools::eps_tolerance<double> tol(get_digits);

            auto G_err = [this](double M) {
//				std::cout << "balanceMach::Calculate()::G_err" << std::endl;
                return balanceG(M,0);
            };

            //std::pair<double, double> res = boost::math::tools::bracket_and_solve_root(G_err, guess, factor, is_rising, tol,
            //                                                                           it);
            std::pair<double, double> res = boost::math::tools::toms748_solve(G_err, Mmin, Mmax, tol, it);

            M2 = res.first + (res.second - res.first) / 2;
            //std::cout << N2_->report() << std::endl;
           // std::cout << "M2: " << M2 << std::endl;
        } else {
           // std::cout << "Using nonLinear" << std::endl;
            M2 = nonLinear<balanceMach>::funSolve(&balanceMach::balanceG, Mmin, Mmax, *this, 0);
            //std::cout << "M2: " << M2 << std::endl;
        }


		double r = balanceG(M2, 0);
		r = r / G1_;
		if (verbosity_ > 0) {
			std::cout << "Normalized Residual = " << r << std::endl;
		}
        return rChok;
	}

//double balanceMach::balanceG(double& M2, const void* par[]) {
//	double dG, dI, dH;
//	delta(dG,dI,dH);
//
//	double G2 = G1_ - dG;
//	double I2 = I1_ - dI;
//	double H2 = H1_ - dH;
//	if(G2==0.0){
//		std::cerr << "Error in Module: " << Name_ << " -- Mass flow at node N2 cannot be zero." << std::endl;
//		exit(1);
//	}
//	if(I2<=0.0){
//		throw std::runtime_error("Error: delta Impulse dI greater than starting impulse I1. Solution impossible.");
//	}
//	double HG2 = H2/G2;
//
//	const void* par1[3];
//	par1[0] = N2_;
//	par1[1] = &M2;
//	par1[2] = &HG2;
//
//	double Tmin = N2_->Tlow();
//	double Tmax = N2_->Thigh()-1500;
//
//	N2_->setPress(N1_->getPress());
//
//    double temp_T = nonLinear<balanceMach>::SolveFalsePos(&balanceMach::THequation,Tmin,Tmax,*this,par1);
//    double temp_P = I2/(1 + N2_->gamma()*std::pow(M2,2.0));
//
//    //N2_->setTP(nonLinear<balanceMach>::SolveFalsePos(&balanceMach::THequation,Tmin,Tmax,*this,par1),I2/(1 + N2_->gamma()*std::pow(M2,2.0)));
//    N2_->setTP(temp_T,temp_P);
//	N2_->setU(M2*N2_->geta());
//
//	return N2_->G() - G2;
//}

//	//double balanceMach::balanceG(double &M2, const void *par[]) {
    double balanceMach::balanceG(double &M2, const void *par[]) {
		double dG, dI, dH;
		delta(dG, dI, dH);

		double G2 = G1_ - dG;
		double I2 = I1_ - dI;
		double H2 = H1_ - dH;
		if (solver_verbosity) {
			fmt::print("\tbalanceMach::balanceG\n");
			fmt::print("\t\tG1_ ({}) - dG ({}) = G2 ({})\n", G1_, dG, G2);
			fmt::print("\t\tI1_ ({}) - dI ({}) = GI ({})\n", I1_, dI, I2);
			fmt::print("\t\tH1_ ({}) - dH ({}) = H2 ({})\n", H1_, dH, H2);
		}
		if (G2 == 0.0) {
			std::cerr << "Error in Module: " << Name_ << " -- Mass flow at node N2 cannot be zero." << std::endl;
			exit(1);
		}
		if (I2 <= 0.0) {
			throw std::runtime_error("Error: delta Impulse dI greater than starting impulse I1. Solution impossible.");
		}
		if (dH > 0.0){

		}
		double HG2 = H2 / G2;

		double guess = 300;
		double factor = 1.25;

		const boost::uintmax_t maxit = 30;
		boost::uintmax_t it = maxit;
		bool is_rising = true;
		int digits = std::numeric_limits<double>::digits;
		int get_digits = digits - 3;
		boost::math::tools::eps_tolerance<double> tol(get_digits);


		N2_->setPress(N1_->getPress());

		auto h_err = [this, M2, HG2](double T) {
			N2_->setTP(T, N2_->getPress());
			if (solver_verbosity) {
				fmt::print("\t\t\tT: {}, p: {}\n", T, N2_->getPress());
				fmt::print("\t\t\t {} - {} = {}\n", N2_->h() + 0.5 * std::pow(M2 * N2_->geta(), 2.0), HG2, N2_->h() + 0.5 * std::pow(M2 * N2_->geta(), 2.0) - HG2);
			}
            return N2_->h() + 0.5 * std::pow(M2 * N2_->geta(), 2.0) - HG2;
		};

		double Tmin = N2_->Tlow();
		double Tmax = N2_->Thigh() - 1500;
		if (solver_verbosity) {
			fmt::print("\t\tTemperature Loop\n");
		}

		std::pair<double, double> r = boost::math::tools::bracket_and_solve_root(h_err, guess, factor, is_rising, tol,
																				 it);

		if (solver_verbosity) {
			if (it >= maxit) {
				fmt::print(
						"\t\tUnable to locate solution in {} iterations: {}\n Current best guess is between {} and {}",
						maxit, it, r.first, r.second);
			} else {
				fmt::print("\t\tConverged after {} iterations (maximum iterations: {})\n", it, maxit);
			}
		}
		double temp_T = r.first + (r.second - r.first) / 2;

		double temp_P = I2 / (1 + N2_->gamma() * std::pow(M2, 2.0));

		//N2_->setTP(nonLinear<balanceMach>::SolveFalsePos(&balanceMach::THequation,Tmin,Tmax,*this,par1),I2/(1 + N2_->gamma()*std::pow(M2,2.0)));

		N2_->setTP(temp_T, temp_P);
        N2_->setU(M2 * N2_->geta());

		if (solver_verbosity) {
			fmt::print("\t\tPressure Loop\n");
		}

		unsigned whileCounter = 0;
		double lastStatus;
        // Check if the Total Enthalpy no longer matches the target Total Enthalpy
		if (fabs(N2_->H() - HG2) > 1e-3 && N2_->equilibrium){
            do {
                // Iterate through the process to find the temperature, with the previous pressure solution
                N2_->setTPX(N1_->getTemp(), temp_P, N1_->X());


                auto h_err = [this, M2, HG2](double T) {
                    N2_->setTP(T, N2_->getPress());
					if (solver_verbosity){
						fmt::print("\t\t\tT: {}, p: {}\n", T, N2_->getPress());
						fmt::print("\t\t\t {} - {} = {}\n", N2_->h() + 0.5 * std::pow(M2 * N2_->geta(), 2.0), HG2, N2_->h() + 0.5 * std::pow(M2 * N2_->geta(), 2.0) - HG2);
					}
                    return N2_->h() + 0.5 * std::pow(M2 * N2_->geta(), 2.0) - HG2;
                };

                std::pair<double, double> r = boost::math::tools::bracket_and_solve_root(h_err, guess, factor, is_rising, tol,
                                                                                         it);

				if (solver_verbosity){
					if (it >= maxit){
						fmt::print("\t\tUnable to locate solution in {} iterations: {}\n Current best guess is between {} and {}", maxit, it, r.first, r.second);
					} else {
						fmt::print("\t\tConverged after {} iterations (maximum iterations: {})\n", it, maxit);
					}
				}
	            temp_T = r.first + (r.second - r.first) / 2;
                temp_P = I2 / (1 + N2_->gamma() * std::pow(M2, 2.0));

				N2_->setTP(temp_T, temp_P);
				N2_->setU(M2 * N2_->geta());
				if (solver_verbosity) {
					fmt::print("\t\tDifference between N2_->H() {}, and HG2 {} is {}\n", N2_->H(), HG2, N2_->H() - HG2);
				}
				if (solver_verbosity) {
					fmt::print("Iterations number: {}, status: {}, progress: {}\n", whileCounter, fabs(N2_->H() - HG2),
							   fabs(lastStatus) - fabs(N2_->H() - HG2));
				}
				lastStatus = N2_->H() - HG2;
				whileCounter++;
            }
            while (fabs(N2_->H() - HG2) > 1e-3);
		}
		if (solver_verbosity) {
			fmt::print("\tDifference between N2_->G() {}, amd G2 {} is {}\n", N2_->G(), G2, N2_->G() - G2);
		}

		return N2_->G() - G2;
	}


//	//double balanceMach::balanceG(double &M2, const void *par[]) {
//	double balanceMach::balanceG(double &M2, const void *par[]) {
//		double dG, dI, dH;
//		delta(dG, dI, dH);
//
//		double G2 = G1_ - dG;
//		double I2 = I1_ - dI;
//		double H2 = H1_ - dH;
//		if (G2 == 0.0) {
//			std::cerr << "Error in Module: " << Name_ << " -- Mass flow at node N2 cannot be zero." << std::endl;
//			exit(1);
//		}
//		if (I2 <= 0.0) {
//			throw std::runtime_error("Error: delta Impulse dI greater than starting impulse I1. Solution impossible.");
//		}
//		if (dH > 0.0){
//
//		}
//		double HG2 = H2 / G2;
//
//		const boost::uintmax_t maxit = 30;
//		boost::uintmax_t it = maxit;
//		int digits = std::numeric_limits<double>::digits;
//		int get_digits = digits - 3;
//		boost::math::tools::eps_tolerance<double> tol(get_digits);
//
//		auto balanceI = [this, M2, HG2, I2](double p){
//			N2_->setTP(N1_->getTemp(), p);
//
//			double guess = 300;
//			double factor = 1.25;
//
//			const boost::uintmax_t maxit = 30;
//			boost::uintmax_t it = maxit;
//			bool is_rising = true;
//			int digits = std::numeric_limits<double>::digits;
//			int get_digits = digits - 3;
//			boost::math::tools::eps_tolerance<double> tol(get_digits);
//
//			auto balanceH = [this, M2, HG2, p](double T){
//				N2_->setTP(T, N2_->getPress());
//				return N2_->h() + 0.5 * std::pow(M2 * N2_->geta(), 2.0) - HG2;
//			};
//
//			std::pair<double, double> r = boost::math::tools::bracket_and_solve_root(balanceH, guess, factor, is_rising, tol,
//																					 it);
//			double temp_T = r.first + (r.second - r.first) / 2;
//			N2_->setTP(temp_T, p);
//			N2_->setU(M2 * N2_->geta());
//
//			return p*(1 + N2_->gamma() * std::pow(M2, 2.0)) - I2;
//		};
//
//		double pmin = 10000;
//		double pmax = 1e8;
//		std::pair<double, double> r = boost::math::tools::toms748_solve(balanceI, pmin, pmax, tol, it);
//
//		double newPressure = r.first + (r.second - r.first) / 2;
//
//		auto balanceH2 = [this, M2, HG2, newPressure ](double T){
//			N2_->setTP(T, newPressure);
//
//			return N2_->h() + 0.5 * std::pow(M2 * N2_->geta(), 2.0) - HG2;
//		};
//
//		double tmax = 6400;
//		double tmin = 300;
//
//		std::pair<double, double> r2 = boost::math::tools::toms748_solve(balanceH2, tmin, tmax, tol, it);
//
//		double newTemp = r2.first + (r2.second - r2.first) / 2;
//		N2_->setTP(newTemp, newPressure );
//		N2_->setU(M2 * N2_->geta());
//
//		fmt::print("T: {} \t P: {}\t U: {}\t M2: {}\n", N2_->getTemp(), N2_->getPress(), N2_->getU(), M2);
//		fmt::print("H: {}\n", N2_->h() + 0.5 * std::pow(M2 * N2_->geta(), 2.0) - HG2);
//		fmt::print("I: {}\n", N2_->getPress()*(1 + N2_->gamma() * std::pow(M2, 2.0)) - I2);
//		fmt::print("G: {}\n", N2_->G() - G2);
//
//
//		return N2_->G() - G2;
//	}


	double balanceMach::THequation(double &T, const void *par[]) {
		Node &N = *(Node *) par[0];
		const double &Mach = *(double *) par[1];
		const double &H = *(double *) par[2];

		//In order to handle instabilities due to the falsePosition algorithm
		if (T <= 0) {
			T = toll_;
		}
		N.setTP(T, N.getPress());
		return N.h() + 0.5 * std::pow(Mach * N.geta(), 2.0) - H;
	}

	double balanceMach::choked() {
		return balanceG(chokingMach_, 0);
	}

	void balanceMach::delta(
			double &dG,
			double &dI,
			double &dH) const {
		dG = 0.0;
		dI = 0.0;
		dH = 0.0;
	}
}
//BOOST_CLASS_EXPORT_IMPLEMENT(balanceMach);

