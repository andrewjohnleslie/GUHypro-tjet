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
 
#ifndef COMPRESSOR_H_
#define COMPRESSOR_H_
 
 
#include "propModule.h"
// include definitions in "propModule.h"

#include "mechanicalLink.h"
#include "Archive.h"
#include "polytropicDuct.h"

 
 namespace hypro {

     class Compressor: public polytropicDuct {

     protected:
 
    	double T04_ ;
    	double T02_;
    	double p02_;
    	mechanicalLink* Shaft_;
 
     public:

    	double pi_c_, effc_, kH_;
		int checkdesignpoint_;


    	Compressor();      	
    	Compressor(std::string name, Node* N1, Node* N2, mechanicalLink* Shaft_, int ID = defID_);
    	Compressor(const Compressor& mod);
    	virtual ~Compressor();

		virtual GPObject &duplicate() { return *(new Compressor(*this)); }
		propModuleSERIAL(Compressor);


     	virtual double calculate();

     	virtual	void   EntropyRise();
		virtual void   PressRise(double &H02);
     	virtual	void   balanceHCompressor();


     	virtual double isoffdesign();
     	double 		   getT04();
     	void           setT04(double& T);
     	void 		   assignForward();

		virtual void   serialize(Archive& ar) const;
	
 
     };


} 

 
 #endif

