
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
 
#include "Turbine.h"
// header file containing variable and function declarations

#include "isentropicDuct.h"
// header that contains isentropic relation (functions) between variables

#include "nonLinear.h"
// header that contains nonLinear solver (function)

#include "Archive.h"

#include <time.h>

namespace hypro {


	Turbine::Turbine()
	:
		polytropicDuct(),
		Shaft_{},
		p02_{0},
		T02_{0},
		efft_{0},
		kH_{0},
		Aratio_{0},
		cp_{0},
		T04_{0},
		checkdesignpoint_{0} {
	}


	Turbine::Turbine(std::string name, Node* N1, Node* N2, mechanicalLink* Shaft_, int ID)
	:
		polytropicDuct(name, N1, N2, ID),
		Shaft_(Shaft_),
		p02_{1.0},
		T02_{1.0},
		efft_{1.0},
		kH_{1.0},
		Aratio_{1.0},
		cp_{1.0},
		T04_{1.0},
		checkdesignpoint_{1} {
	}


	Turbine::Turbine(const Turbine& mod)
	:
		polytropicDuct(mod),
		Shaft_(mod.Shaft_),
		p02_(mod.p02_),
		T02_(mod.T02_),
		efft_(mod.efft_),
		kH_(mod.kH_),
		Aratio_(mod.Aratio_),
		T04_(mod.T04_),
		cp_(mod.cp_),
		checkdesignpoint_(mod.checkdesignpoint_) {
	}

	Turbine::~Turbine() {
	}
	

	double Turbine::calculate() {
	cnt1++;
	clock_t start1 = clock();
		double rChok = 1.0;

		if (Xswitch_ == 0){ 	// Frozen Flow conditions,
			N2_->X(N1_->X());
		}
		if (Xswitch_ == 1){
			N2_->X(XSET_);
		}
		balanceHTurbine();
		assignForward();
	clock_t end1 = clock();
	elapsed1 = double(end1 - start1)/CLOCKS_PER_SEC;
	
	//std::cout << cnt1 << ": CPU time measured (Turbine::calculate): " << elapsed1 << " s"<< std::endl;	
		
		return rChok;

	}


	void Turbine::balanceHTurbine() {
        cnt2++;
	clock_t start2 = clock();
		double p01 = N1_->getp0();
		double T01 = N1_->getT0();

		double cpmin = 1000;  //Taken from cp(gasmodel.Tlow() for Jet A exhaust (actual value is ~1032)
		//double cpmax = 1220;
		double cpmax = N1_->Cp(); //The value for Cp will be lower at the exit due to the expansion process

		const void *parbalancecp[2] = {&T01, &p01};
		nonLinear<Turbine>::SolveFalsePos(&Turbine::balanceCpTurbine, cpmin, cpmax, *this, parbalancecp);
	clock_t end2 = clock();
	elapsed2 = double(end2 - start2)/CLOCKS_PER_SEC;
	
	//std::cout << cnt2 << ": CPU time measured (Turbine::balanceHTurbine): " << elapsed2 << " s"<< std::endl;
	}



	double Turbine::balanceCpTurbine(double& cp, const void *par[]) {
        cnt3++;
	clock_t start3 = clock();
		double T01 = *(double *) (par[0]);
		double p01 = *(double *) (par[1]);	

		double Power = Shaft_->getPower();	
		double R = N1_->R()/N1_->W();
		double H2s;
		Node Nt(*N2_);
		if (this->checkdesignpoint_ == 1.0) {
//			******************* For .H() method ***************************************************************
			double deltaH = Power/N1_->mfr();	
			H2s = N1_->H()  - (deltaH);	
			Nt.setHP(H2s,p01);
			Nt.setU(0);
			T02_ = Nt.getTemp();
//			***************************************************************************************************
		}

		else{

			kH_ = getkH(cp);
			T02_ = (T01 * (1 - kH_));
			Nt.setTP(T02_,p01);
			Nt.setU(0);
//			Nt.printOut(std::cout,true);
			H2s = Nt.H();		
		}

		PressDrop(H2s);

		double r = (cp - N2_->Cp());
		std::cout << "cp1 is: " << cp << "\n";
		
		clock_t end3 = clock();
		elapsed3 = double(end3 - start3)/CLOCKS_PER_SEC;
	
		//std::cout << cnt3 << ": CPU time measured (Turbine::balanceHTurbine): " << elapsed3 << " s"<< std::endl;
		return r;

	}

	




	void Turbine::PressDrop(double &H02_){
	cnt4++;
	clock_t start4 = clock();
			double dp0, dH0, H0x;

			double R = N1_->R()/N1_->W();
			Node Nt(*N2_);
			Nt.setU(0);
			Nt.A(N2_->A());

	        double T0x = N1_->getT0();
			double p0x = N1_->getp0();
			H0x = N1_->H();
	        dH0 = (H0x - H02_)/1000;		//ORG: 100000

			for (int i = 0; i < 1000; i++){	//ORG: 100000

				dp0 = p0x * (dH0/(T0x*R*efft_));
				p0x -= dp0;
				H0x -= dH0;
				T0x = Nt.getTemp();
				Nt.setHP(H0x,p0x);

			}

			p02_ = p0x;
			Nt.setHP(H0x,p0x);
			double s0x = Nt.s();

//			std::cout<< "H0x is:" << H0x <<'\n';
//			std::cout<< "T0x is:" << T0x <<'\n';
//			std::cout<< "p0x is:" << p0x <<'\n';

//			if (this->checkdesignpoint_ == 1.0) {
//			Nt.printOut(std::cout, true);
//			}

			double Tmin {800};
			double Tmax {T0x};
			const void *passHP[] = {&s0x, &H0x};
			nonLinear<Turbine>::SolveFalsePos(&Turbine::HmodeT, Tmin, Tmax, *this, passHP);

		clock_t end4 = clock();
		elapsed4 = double(end4 - start4)/CLOCKS_PER_SEC;
	
		//std::cout << cnt3 << ": CPU time measured (Turbine::balanceHTurbine): " << elapsed4 << " s"<< std::endl;

		}


	void Turbine::assignForward(){
		double T01 = N1_->getT0();

		if (this->checkdesignpoint_) {
			kH_ =   1 - (T02_/T01);
			double exponent = (((2*efft_)*(N1_->gamma()-1))/((2*N1_->gamma())-((efft_)*(N1_->gamma()-1))));
			double Arcp = std::pow((1 - kH_),(1/exponent));
			Aratio_ = Arcp*(std::pow((N1_->Cp()/N2_->Cp()),0.5));
			setT04(T01);
			std::cout << "kH is: " << kH_ <<"\n";
			
		}
		else {
			setT04(T01);
			Shaft_->mfr_ = N1_->mfr();
			Shaft_->Cp2_ = N2_->Cp();
			Shaft_->Cp_  = N1_->Cp();
			Shaft_->Ht_  = N1_->mfr()*(N1_->H() - N2_->H());
			Shaft_->kH_  = kH_;


		}
	}


	double Turbine::isoffdesign() {
		return checkdesignpoint_;
	}


	double Turbine::getT04(){
		return T04_;
	}


	void Turbine::setT04(double& T){
		T04_ = T;
	}


	double Turbine::getkH(double &cp){
		double exponent = (((2*efft_)*(N1_->gamma()-1))/((2*N1_->gamma())-((efft_)*(N1_->gamma()-1))));
		double cpratio = std::pow((cp/N1_->Cp()),0.5);
		std::cout << "cp is: " << cp << "\n";
		kH_ = 1 - (std::pow((Aratio_*cpratio),exponent));
		//std::cout << "kH in " << N2_->Name_ << " is: " << kH_ <<"\n";			// UNCOMMENT
		//std::cout << "Aratio_ in " << N2_->Name_ << " is: " << Aratio_ <<"\n";		// UNCOMMENT	
		return kH_;
	}


	void Turbine::serialize(Archive& ar) const{
		propModule::serialize(ar);
		
		ar.put("Area Ratio", Aratio_);
		ar.put("Off Design Constant kH", kH_);
		ar.put("Turbine Efficiency", efft_);
		ar.put("Mechanical Link ID",Shaft_->ID_);

	}

}
	
	
