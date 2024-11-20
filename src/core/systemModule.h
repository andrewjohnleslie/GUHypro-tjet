/*!
 *  \brief      Module representing the whole propulsion system.
 *  \details    It contains all the system components.
 *              Execute all the required iterations to solve the whole model.
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

#ifndef SYSTEMMODULE_H_
#define SYSTEMMODULE_H_

#include "Collection.h"

namespace hypro {
	class systemModule : public Collection {
		///Function that is solved in case of feedback.
		/**@param[in] x value of the reduce parameter of the feedback destination
         * @param[in] par array of parameters.
         * 	- `par[0]` is an iterator pointing to the feedback destination.
         * 	- `par[1]` is an iterator pointing to the feedback source.
         * @return residual of the feedback source calculation
         *
         * The function first set the reduce parameter and then
         * it evaluates in sequence all the modules between
         * the feedback destination and its origin.
         * It eventually returns the residual outcoming from the feedback
         * origin.
         */
		double findChok     (double &x, const void *par[]);
		double offError_;
		const double offToll_ = 0.05;
		double Tmin_ {1000};
		double Tmax_ {1700};

	public:

		typedef std::list<std::pair<ModulePtr,ModulePtr>> OffDesignPairedList;
		OffDesignPairedList PairedList_;

		///Changes the behaviour of the solver
		/**If true the feedback loop search limits are reset every time
         * a choking is detected. This slows down the calculation but might be more robust.
         * This option can in same cases generates choking in the choking find loop, throwing an error.
         * The option is false by default.
         */
		static bool resetLoopLimits_;

		///Constructor with name
		/**@param[in] N1 input node
         * @param[in] N2 output node
         * @param[in] name name of the module
         *
         * The name of the module will be initialised to the default one.
         */
		systemModule(std::string name, Node *N1, Node *N2);

		systemModule(const systemModule &mod);

		systemModule() {}

		virtual ~systemModule(); ///<Destructor

		virtual GPObject &duplicate() { return *(new systemModule(*this)); }

		//propModuleSERIAL(systemModule);
		using Collection::serialize;

		virtual double calculate();

		///Handles the error solving and iteration of off design characteristics
        double offdesigncalculate(const void *par[]);

		///Checks and identifies any modules running off design
		std::shared_ptr<propModule> offdesigncheck (const void *par[]);

		///Handles the error solving and iteration of off design characteristics
        void offdesignlistcalculate(OffDesignPairedList PairedList);

		///Function used for root finding method of Off Design Temperature
        double offdesignTemperatureSeek(double &T, const void *par[]);

		///Returns the total thrust of the propulsion system `[N]`
		double thrust() const;

		///Prints the thrust that any subtree of the system would generate if connected with an ideal adapted nozzle.
		void thrustPath();

		Glib::RefPtr<ModuleGraph> draw();

		///Calculate the total cross sectional area of the engine.
		/***The area is defined as the sum of the inlet area of
         * all the intakes of the engine.
         */
		double crossArea() const;

		void save(std::ofstream &fs);

		static systemModule *load(std::ifstream &fs);

		static void deleteAll(systemModule *);

        virtual void serialize(Archive& ar)const;
	};
}
//BOOST_CLASS_EXPORT_KEY(systemModule);

#endif /* SYSTEMMODULE_H_ */
