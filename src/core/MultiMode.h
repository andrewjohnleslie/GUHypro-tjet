/*!
 *  \brief      System model for multi mode combined cycle engines
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

#ifndef MULTIMODE_H_
#define MULTIMODE_H_

#include "propModule.h"
#include "systemModule.h"

namespace hypro {
	class MultiMode : public propModule {

	protected:
		int presentMode_;    ///<present selected mode
		int preSwitchMode_;   ///<previously selected mode

		///Defines the priority for mode selection
		/**When two or more modes could be operating at a certain condition.
         * This priority allows a choice.
         *
         * @return the list of mode indexes ordered from high to low priority
         */
		std::list<int> priority() const;

		MultiMode() {}

	public:

		std::vector<systemModule> modes_; ///<Vector containing the models of each propulsive mode

		MultiMode(std::string name, Node *N1, Node *N2);

		MultiMode(const MultiMode &mod);

		~MultiMode();

		virtual double calculate();

		using propModule::N1;

		//Make sure that whenever N1_ is changed the N1_ of all the modes is updated
		void N1(Node &N1);

		using propModule::N2;

		//Make sure that whenever N2_ is changed the N2_ of all the modes is updated
		void N2(Node &N2);

		///Automatically select the propulsion mode
		/**The function select the propulsion mode based typically on free stream conditions.
         * It has to be implemented in subclasses
         *
         * @return index of the selected mode
         */
		virtual int autoSelect() const =0;

		///Switch event function for Matlab ode solver
		/**Event function useful for inclusion within trajectory integration.
         *
         * @return negative number if the switching threshold has been overtaken
         */
		virtual double events() const =0;

		///Allows the module to switch mode during the calculation
		void switchMode();

		///Adds a new empty mode
		/**
         * @param Name name of the mode to be added
         * @return a reference to the mode just added
         */
		virtual systemModule &add(std::string Name);

		///Adds a new mode copy of an existing mode
		/**
         * @param Name name of the new mode
         * @param i index of the existing mode that will be copied
         * @return a reference to the mode just added
         */
		virtual systemModule &add(std::string Name, unsigned i);

		///Calculate the thrust of the engine
		double thrust() const;

		std::vector<double> propellantMfr() const;

		void unreduce();
	};
}
#endif /* MULTIMODE_H_ */
