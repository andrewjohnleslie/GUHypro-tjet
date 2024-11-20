
/*!
// Created by Mark De Luca on 07/12/18.
//
// For any queries, contact m.de-luca.1@research.gla.ac.uk
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
 * along with HyPro.  If not, see <http://www.gnu.org/licenses/>. 
*/
 
 
#ifndef TURBINEBLEED_H_
#define TURBINEBLEED_H_
 

#include "Turbine.h"
// include definitions in "propModule.h"

#include "mechanicalLink.h"
#include "Archive.h"

namespace hypro {
 
 	 class TurbineBleed: public Turbine {

 	 protected:

 		 Node *N3_;

	 public:


		TurbineBleed();

		TurbineBleed(std::string name, Node* N1, Node* N2, Node* N3, mechanicalLink* Shaft_, int ID = defID_);

		TurbineBleed(const TurbineBleed& mod);

		virtual ~TurbineBleed();

		virtual GPObject &duplicate() { return *(new TurbineBleed(*this)); }


		virtual double 	calculate();
		virtual void 	changeComposition();
		virtual void   	ondesign();
		virtual double 	balanceCpTurbine(double& T, const void* par[]);
		virtual double 	balanceGTurbine(double& M2, const void* par[]);
		virtual double 	balanceHBleed(double& M2, const void* par[]);

 
    };
  
 
 
 }
 
 #endif
