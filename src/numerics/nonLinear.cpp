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

#include "nonLinear.h"
#include <iostream>
#include <cmath>
#include <stdexcept>

namespace hypro {
	class isentropicDuct;

	class balanceMach;

	class systemModule;

	class Wedge;

	class AdaptedNozzle;

	class MyGene;

	class thermoState;

	class thermoKineticState;

	class Fanno;

    class Compressor;

    class Turbine;

    class TurbineBleed;

    class polytropicDuct;


	template<class Module>
	const double nonLinear<Module>::toll_ = 1e-6;

	template<class Module>
	const double nonLinear<Module>::hightoll_ = 1e-12	;


	template<class Module>
	int nonLinear<Module>::maxIter_ = 10000;

	template<class Module>
	nonLinear<Module>::nonLinear() {
		// TODO Auto-generated constructor stub

	}

	template<class Module>
	nonLinear<Module>::~nonLinear() {
		// TODO Auto-generated destructor stub
	}

	template<class Module>
	double nonLinear<Module>::funSolve(double(Module::*func)(double &x, const void *par[]),
									   double x1, double x2, Module &obj, const void *par[]) {


		double y1 = (obj.*func)(x1, par);
		double y2 = (obj.*func)(x2, par);

		if (y1 * y2 > 0) {
			throw std::runtime_error("Error: function evaluations at starting points must be of opposite signs.");
		}

		double r = 2 * toll_;
		int i = 0;
		double x;
		while (r > toll_ and i < maxIter_) {
			x = 0.5 * (x1 + x2);

			double y = (obj.*func)(x, par);

			if (y * y1 < 0) {
				y2 = y;
				x2 = x;
			} else {
				y1 = y;
				x1 = x;
			}

			r = std::abs(x1 - x2);
			i++;
		}

		if (i == maxIter_) {
			std::cout << "Warning: maximum number of iteration reached." << std::endl;
		}

		return x;
	}

	template<class Module>
	double nonLinear<Module>::funSolveHighTol(double(Module::*func)(double &x, const void *par[]),
										   double x1, double x2, Module &obj, const void *par[]) {


			double y1 = (obj.*func)(x1, par);
			double y2 = (obj.*func)(x2, par);

			if (y1 * y2 > 0) {
				throw std::runtime_error("Error: function evaluations at starting points must be of opposite signs.");
			}

			double r = 2 * hightoll_;
			int i = 0;
			double x;
			while (r > hightoll_ and i < maxIter_) {
				x = 0.5 * (x1 + x2);

				double y = (obj.*func)(x, par);

				if (y * y1 < 0) {
					y2 = y;
					x2 = x;
				} else {
					y1 = y;
					x1 = x;
				}

				r = std::abs(x1 - x2);
				i++;
			}

			if (i == maxIter_) {
				std::cout << "Warning: maximum number of iteration reached." << std::endl;
			}

			return x;
		}



	template<class Module>
	double nonLinear<Module>::SolveFalsePos(double(Module::*func)(double &x, const void *par[]),
											double x1, double x2, Module &obj, const void *par[]) {

		double y1 = (obj.*func)(x1, par);
		double y2 = (obj.*func)(x2, par);

		//std::cout << y1 << std::endl;
		//std::cout << y2 << std::endl;

		if (y1 == y2) {
			std::cout << "Error" << std::endl;
			throw std::runtime_error("Error: function evaluations at starting points must be different.");
		}


		double r = 2 * toll_;
		int i = 0;
		double x;
		while (r > toll_ and i < maxIter_) {
			x = x1 + y1 * (x1 - x2) / (y2 - y1);

			double y = (obj.*func)(x, par);

			r = std::min(std::abs(x - x2), std::abs(x - x1));

			if (y * y1 < 0) {
				y2 = y;
				x2 = x;
			} else {
				y1 = y;
				x1 = x;
			}
			i++;
		}

		if (i == maxIter_) {
			std::cout << "Warning: maximum number of iteration reached." << std::endl;
		}
		//std::cout << "NL::SFP " << "x: " << x << std::endl;
		return x;
	}

	template<class Module>
	void nonLinear<Module>::printFunc(double(Module::*func)(double &x, const void *par[]),
									  double x1, double x2, Module &obj, const void *par[], int N) {

		double dx = (x2 - x1) / N;
		for (double x = x1; x <= x2; x += dx) {
			std::cout << x << " " << (obj.*func)(x, par) << std::endl;
		}
	}

//Template implementation
	template
	class nonLinear<isentropicDuct>;

	template
	class nonLinear<balanceMach>;

	template
	class nonLinear<systemModule>;

	template
	class nonLinear<Wedge>;

	template
	class nonLinear<AdaptedNozzle>;

	template
	class nonLinear<MyGene>;

	template
	class nonLinear<thermoState>;

	template
	class nonLinear<thermoKineticState>;

	template
	class nonLinear<Fanno>;

	template
	class nonLinear<Compressor>;

	template
	class nonLinear<Turbine>;

	template
	class nonLinear<TurbineBleed>;

	template
	class nonLinear<polytropicDuct>;
}
