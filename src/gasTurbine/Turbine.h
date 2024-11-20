
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
 
 
#ifndef TURBINE_H_
#define TURBINE_H_
 

#include "propModule.h"
// include definitions in "propModule.h"

#include "mechanicalLink.h"
#include "Archive.h"
#include "polytropicDuct.h"

namespace hypro {
 
 	 class Turbine: public polytropicDuct {
 
 	 protected:
		int cnt1 = 0;
		int cnt2 = 0;
		int cnt3 = 0;
		int cnt4 = 0;
  		double cp_;
 		double T04_;
		double T02_;
		double p02_;
		mechanicalLink* Shaft_;
		double elapsed1 = 0;
		double elapsed2 = 0;
		double elapsed3 = 0;
		double elapsed4 = 0;

	 public:

		double kH_;
		double Aratio_;
		double efft_;
		int checkdesignpoint_;
		int Xswitch_ {0};
		std::vector<double> XSET_;





		Turbine();

		Turbine(std::string name, Node* N1, Node* N2, mechanicalLink* Shaft_, int ID = defID_);

		Turbine(const Turbine& mod);

		virtual ~Turbine();

		virtual GPObject &duplicate() { return *(new Turbine(*this)); }

		propModuleSERIAL(Turbine);

		virtual double calculate();
		virtual void   balanceHTurbine();
		virtual double balanceCpTurbine(double& T, const void* par[]);
		virtual void   PressDrop(double &H02);


		virtual void   assignForward();
		virtual double isoffdesign();
     	double 		   getT04();
     	void           setT04(double& T);
     	virtual double getkH(double &cp);


		virtual void   serialize(Archive& ar) const;


 
 
    };
  
 
 
 }
 
 #endif
