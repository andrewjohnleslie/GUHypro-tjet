/*!
 *  \brief      Provides methods to solve non linear root searching problems
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

#ifndef NONLINEAR_H_
#define NONLINEAR_H_

namespace hypro {
	template<class Module>
	class nonLinear {

		static const double toll_; ///<Numeric tolerance;
		static const double hightoll_; ///<Numeric tolerance;


	public:
		static int maxIter_; ///<Maximum number of iterations

		nonLinear();

		virtual ~nonLinear();

		///Bisection method
		/**
         * @param func pointer to the function to be solved
         * @param x1 lower boundary for the searched variable
         * @param x2 upper boundary for the searched variable
         * @param obj reference to the object for which the method func is called
         * @param par parameter vector
         * @return
         */

		static double funSolve(double(Module::*func)(double &x, const void *par[]),
							   double x1, double x2, Module &obj, const void *par[]);
		///Bisection method
				/**
		         * @param func pointer to the function to be solved
		         * @param x1 lower boundary for the searched variable
		         * @param x2 upper boundary for the searched variable
		         * @param obj reference to the object for which the method func is called
		         * @param par parameter vector
		         * @return
		         */
		static double funSolveHighTol(double(Module::*func)(double &x, const void *par[]),
									   double x1, double x2, Module &obj, const void *par[]);

		///Secant method with False Position
		/**
         * @param func pointer to the function to be solved
         * @param x1 lower boundary for the searched variable
         * @param x2 upper boundary for the searched variable
         * @param obj reference to the object for which the method func is called
         * @param par parameter vector
         * @return
         */
		static double SolveFalsePos(double(Module::*func)(double &x, const void *par[]),
									double x1, double x2, Module &obj, const void *par[]);

		///Print function at an uniform grid of N points
		static void printFunc(double(Module::*func)(double &x, const void *par[]),
							  double x1, double x2, Module &obj, const void *par[], int N);
	};
}
#endif /* NONLINEAR_H_ */
