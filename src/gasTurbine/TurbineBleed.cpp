
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
 
#include "TurbineBleed.h"
// header file containing variable and function declarations

#include "isentropicDuct.h"
// header that contains isentropic relation (functions) between variables

#include "nonLinear.h"
// header that contains nonLinear solver (function)

#include "Archive.h"

namespace hypro {


	TurbineBleed::TurbineBleed()
	:
		Turbine(),
		N3_(NULL){
	}


	TurbineBleed::TurbineBleed(std::string name, Node* N1, Node* N2, Node* N3, mechanicalLink* Shaft_, int ID)
	:
		Turbine(name, N1, N2, Shaft_, ID),
		N3_(N3){
	}


	TurbineBleed::TurbineBleed(const TurbineBleed& mod)
	:
		Turbine(mod),
		N3_(mod.N3_) {
	}


	TurbineBleed::~TurbineBleed() {
	}
	
	

	double TurbineBleed::calculate() {

		changeComposition();
		ondesign();


		return 1.0;

	}

	void TurbineBleed::changeComposition() {
		double N2 = N1_->Nfr() + N3_->Nfr();
		N2_->X((N1_->Nfr() * N1_->X() + N3_->Nfr() * N3_->X()) / N2);
	}


	void TurbineBleed::ondesign() {


		double p01 = N1_->getp0();
		double T01 = N1_->getT0();

		double cpmin = 1000;  //Taken from cp(gasmodel.Tlow() for Jet A exhaust (actual value is ~1032)
		double cpmax = N1_->Cp()/1.2; //The value for Cp will be lower at the exit due to the expansion process


		const void *parbalancecp[2] = {&T01, &p01};
		double cp2 = nonLinear<TurbineBleed>::funSolve(&TurbineBleed::balanceCpTurbine, cpmin, cpmax, *this, parbalancecp);



//		kH_ =   1 - (T02_/T01);
//		double exponent = (((2*efft_*N1_->gamma())*(N1_->gamma()-1))/((2*N1_->gamma())-((efft_)*(N1_->gamma()-1))));
//		Aratio_ = std::pow(kH_,(1/exponent));
//		setT04(T01);
//
//		std::cout<<N1_->s()<<std::endl;
//		std::cout<<N2_->s()<<std::endl;


	}
	

	double TurbineBleed::balanceCpTurbine(double& cp, const void *par[]) {


//			std::cout << "In Function: " << BOOST_CURRENT_FUNCTION << std::endl;

			double T01 = *(double *) (par[0]);


			double Ptmin = 1.3*N1_->getp0();
			double sf = (N1_->mfr() - 11)*3;
			double Ptmax = Ptmin/sf;
			const void *par2[2] = {&T01, &cp};


			double Pt = nonLinear<TurbineBleed>::funSolve(&TurbineBleed::balanceHBleed, Ptmin, Ptmax, *this, par2);


			double r = (cp - N2_->Cp());
			return r;




	}


	double TurbineBleed::balanceHBleed(double& Pt, const void *par[]) {


		// First, guess the power extract for the bleed portion and use to calculate p02
		//////////////////////////////////////////////////////


//		std::cout << "In Function: " << BOOST_CURRENT_FUNCTION << std::endl;

		double T01 = *(double *) (par[0]);
		double p01 = N1_->getp0();
		double cp =  *(double *) (par[1]);

		double Power = Shaft_->getPower();
		Node N1x(*N1_);
		Node N3x(*N3_);

        p02_ = Pt;
		N1x.isentropicP(*N1_, Pt);
		N3x.isentropicP(*N3_, Pt);

	    double ht_out_c = ((N1_->mfr()/(N1_->mfr()+N3_->mfr())) * ((N1_->getT0() * N1_->Cp()) + (efft_*(N1x.getT0() * N1x.Cp() - N1_->getT0() * N1_->Cp())))); 
        double ht_out_b = ((N3_->mfr()/(N1_->mfr()+N3_->mfr())) * ((N3_->getT0() * N3_->Cp()) + (efft_*(N3x.getT0() * N3x.Cp() - N3_->getT0() * N3_->Cp()))));

        double ht_out = ht_out_c + ht_out_b;
        
        T02_ = ht_out/cp;





		double Mmin = 0.0;
		double Mmax = 1.0;	
		const void *par2[2] = {&p02_, &T02_};

		double M2 = nonLinear<TurbineBleed>::funSolve(&TurbineBleed::balanceGTurbine, Mmin, Mmax, *this, par2);


		double PowerReq = ((efft_*(N1_->getT0() * N1_->Cp() - N1x.getT0() * N1x.Cp())) * N1_->mfr()) + ((efft_*(N3_->getT0() * N3_->Cp() - N3x.getT0() * N3x.Cp())) * N3_->mfr());
		double r = N2_->getp0() - Pt;
		return r;

	}



	double TurbineBleed::balanceGTurbine(double& M2, const void *par[]) {

//		std::cout << "In Function: " << BOOST_CURRENT_FUNCTION << std::endl;


		double p2_p02, T2_T02;
		isentropicDuct::isentropic(M2, N1_->gamma(), p2_p02, T2_T02);

		double T2, U2, p2;

		T2 = (T02_ * T2_T02);
		N2_->setTemp(T2);

		U2 = (M2 * N2_->geta());
		N2_->setU(U2);

		p2 = (p02_ * p2_p02);
		N2_->setPress(p2);

		double r = (N2_->mfr() - (N1_->mfr() + N3_->mfr()));
	    return r;


	}

/*
// Calculate the enthalpic ideal temperature from mixing and extracting, then calculate pressure with polytropic efficiency to work out actual pressure;
		T02_ = ((N1_->mfr() * T01 * N1_->Cp()) + (N3_->mfr() * N3_->getT0() * N3_->Cp()) - (Power + BldEnth)) / (cp * (N1_->mfr() + N3_->mfr()));
		double tau_t0 = T02_/T01;
		p02_ = std::pow(tau_t0,(N1_->gamma())/((N1_->gamma() - 1 ) * efft_)) * N1_->getp0();

		// calculate T02 and thus p02 as if no cooling, then to be later modified
		double T02x =  ((T01*N1_->Cp()) - (Power / N1_->mfr()))/cp;
		double T02x_T01 = T02x/T01;
		double p02x = std::pow(T02x_T01,(N1_->gamma())/((N1_->gamma() - 1 )* efft_ )) * N1_->getp0();
		p02_= p02x;

		// Next, find the value for p2 for the calculated  value of p02, and then calculate the ennthalpy of this flow from the bleed by isentropicly expanding bleed flow to this pressure
		//////////////////////////////////////////////////////
*/

}
	
	
