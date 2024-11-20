
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

#include "Compressor.h"
// header file containing variable and function declarations

#include "isentropicDuct.h"
// header that contains isentropic relation (functions) between variables

#include "nonLinear.h"
// header that contains nonLinear solver (function)

#include "Archive.h"
#include <iostream>

namespace hypro {


    Compressor::Compressor()
    :
        polytropicDuct(),
		Shaft_{},
		p02_{0},
		T02_{0},
		T04_{0},
		pi_c_{0},
		effc_{0},
		kH_{0},
		checkdesignpoint_{0} {
    }
    
    Compressor::Compressor(std::string name, Node* N1, Node* N2, mechanicalLink* Shaft_, int ID)
	:
		polytropicDuct(name, N1, N2, ID),
		Shaft_(Shaft_),
		p02_{1.0},
		T02_{1.0},
		T04_{1.0},
		pi_c_{1.0},
		effc_{1.0},
		kH_{1.0},
		checkdesignpoint_{1} {
	}

	Compressor::Compressor(const Compressor& mod)
	:
		polytropicDuct(mod),
		Shaft_(mod.Shaft_),
		p02_(mod.p02_),
		T02_(mod.T02_),
		T04_{mod.T04_},
		pi_c_(mod.pi_c_),
		effc_(mod.effc_),
		kH_(mod.kH_),
		checkdesignpoint_(mod.checkdesignpoint_) {
	}

	Compressor::~Compressor(){
	}
	

	double Compressor::calculate(){

		double rChok = 1.0;
		N2_->X(N1_->X());
		
		balanceHCompressor();
		assignForward();
		
		return rChok;

	}


	void Compressor::balanceHCompressor(){


		double T01 = N1_->getT0();
		double p01 = N1_->getp0();

		if(this->checkdesignpoint_ == 1.0){
			p02_ = N1_->getp0() * pi_c_;		// pi_c usage - design
			EntropyRise();
		}

		else{
			kH_ = Shaft_->kH_;
			Node Nc(*N2_);
			double deltaH = Shaft_->Ht_/N1_->mfr();	
			double H02 = N1_->H() + deltaH;		
			Nc.setU(0);
			Nc.setHP(H02,p01);		
			T02_ = Nc.getTemp();
//			std::cout << "\nT02 is: " << T02_ << "\n";
//			Nc.printOut(std::cout, true);
			PressRise(H02);
		}
		return;

	}



	void Compressor::EntropyRise(){

		Node Nt(*N1_);
		double dp0 = (N1_->getp0() * (pi_c_ - 1)) /1000;	//ORG: 100000		// pi_c usage - design
		double ds0 {0};

		double s0x = N1_->s();
		double p0x = N1_->getp0();
		double R = N1_->R()/N1_->W();


		for (int i = 0; i < 1000; i++){ 	//ORG: 100000

			p0x = p0x + dp0;
			ds0 = R * ((1-effc_)/effc_) * (dp0/p0x);
			s0x += ds0;
			
		}

		// Create a dummy node to solve for calculate entropy condition if p02 was static pressure. This defines the upper limit of T2 for the PmodeC numerical solve.
		Nt.setSP(s0x,p0x);

		double Tmin {400};
		double Tmax {Nt.getTemp()};
		const void *passSP[] = {&s0x, &p0x};
//		std::cout<< "R is:" << R <<'\n';

		nonLinear<Compressor>::SolveFalsePos(&Compressor::PmodeC, Tmin, Tmax, *this, passSP);

	}



	void Compressor::PressRise(double &H02_){

		double dp0, dH0;
		double R = N1_->R()/N1_->W();
		Node Nt(*N1_);
		Nt.setU(0);
		Nt.A(N1_->A());
//		double ds0 {0};

//		double s0x = N1_->s();
		double T0x = N1_->getT0();
		double p0x = N1_->getp0();
		double H0x = N1_->H();
        dH0 = (H02_ - H0x)/10000;	//ORG: 1000000
		

		for (int i = 0; i < 10000; i++){	//ORG: 1000000

			dp0 = p0x * (effc_) * (dH0/(T0x*R));
			H0x += dH0;
			p0x = p0x + dp0;
			T0x = Nt.getTemp();
			Nt.setHP(H0x,p0x);
			Nt.setU(0);

		}

		dp0 = p0x * (effc_) * (dH0/(T0x*R));
		p0x = p0x + dp0;

		p02_ = p0x;
		Nt.setHP(H0x,p0x);
		double s0x = Nt.s();

//		std::cout<< "H0x is:" << H0x <<'\n';		
//		std::cout<< "T0x is:" << T0x <<'\n';
//		std::cout<< "p0x is:" << p0x <<'\n';		
//		std::cout<< "s0x is:" << s0x <<'\n';

//		Nt.printOut(std::cout, true);

		double Tmin {400};
		double Tmax {T0x};
		const void *passHP[] = {&s0x, &p0x};

		nonLinear<Compressor>::SolveFalsePos(&Compressor::PmodeC, Tmin, Tmax, *this, passHP);

	}





	void Compressor::assignForward(){

		if (this->checkdesignpoint_) {
			double EnthalpyDifference =   N1_->mfr()*( N2_->H() - N1_->H());

			Shaft_->setPower(EnthalpyDifference);
			Shaft_->mfr_ = N1_->mfr();
		}
		else {

			pi_c_ = N2_->getp0()/N1_->getp0(); //This line updates the member to store the new calculated value, for use when serialized @ off design conditions

		}
	}


	double Compressor::isoffdesign() {
		return checkdesignpoint_;
	}

	double Compressor::getT04(){
		return T04_;
	}

	void Compressor::setT04(double& T){
		T04_ = T;
	}

	void Compressor::serialize(Archive& ar) const{
		propModule::serialize(ar);
		
		ar.put("Compressor Ratio", pi_c_);
		ar.put("Compressor Efficiency", effc_);
		ar.put("Mechanical Link ID",Shaft_->ID_);

	}

}
