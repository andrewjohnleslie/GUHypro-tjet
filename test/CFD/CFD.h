/*!
 *  \brief      Validation test against CFD simulation of the Marquardt SERJ engine
 *  \details    See section 4.5 in <a href="http://oleg.lib.strath.ac.uk/R/?func=dbin-jump-full&object_id=26892">
 *              "Toward automated design of Combined Cycle Propulsion." by Alessandro Mogavero</a>
 *
 *              @param selCase Test case identification
 *              @param eta mixer efficiency
 *              @return pointer to the system module
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

#pragma once

#define DO_QUOTE(X)       #X
#define QUOTE(X)          DO_QUOTE(X)

#include <iostream>
#include <time.h>
#include <map>
#include <fenv.h>

#include "core/Node.h"

//#include <IOdictionary.H>
//#include <IOobject.H>
//#include <Time.H>
//#include <argList.H>
//#include <specieThermo.H>
//#include <janafThermo.H>
//#include <perfectGas.H>

//using namespace Foam;

//#include <mixture.h>
#include "solvers/isentropicDuct.h"
#include "core/systemModule.h"
#include "combustion/Combustion.h"
#include "injectors/InjectionPhi.h"
#include <AdaptedThroatInlet.h>
#include "Wedge.h"
#include "Rocket.h"
#include "Mixer.h"
#include <AdaptedNozzle.h>
#include "core/MultiModeMachP0.h"
#include <AdaptedInlet.h>
#include "NeutralLink.h"
#include <chokedConv.h>

#include <boost/math/tools/roots.hpp>
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/Phase.h"
#include <memory>

// using namespace Cantera;

enum caseTag{
	CASE1 = 0,
	CASE1tuned = 1,
	CASE2 = 2
};

std::unique_ptr<hypro::systemModule> CFD(caseTag selCase, double eta){

	double foot =  0.3048; //conversion feet to meters
	double lbm = 0.45359; //conversion libre to kg
	double psi = 6894.75728; //conversion psi to Pa
	double btu = 1055.06; //conversion btu to J
	double pi = 3.14159265359;

	// Thermal model
	//Time runTime(QUOTE(HyProRescource),"");

	//IOdictionary a
	//(
	//		IOobject
	//		(
	//				"thermophysicalProperties",
	//				"",
	//				runTime,
	//				IOobject::MUST_READ,
	//				IOobject::NO_WRITE,
	//				false
	//		)
	//);
	//wordList species(4,"H2");
	//species[1] = "O2";
	//species[2] = "N2";
	//species[3] = "H2O";

	//mixture< specieThermo< janafThermo<perfectGas> > >& specieGas =
	//		*new mixture< specieThermo< janafThermo<perfectGas> > >(a,species);

	//List<double> X(4,0.0);
	//X[1] = 0.209;   //Air
	//X[2] = 1-0.209;

	hypro::Node freeStream("/usr/local/share/cantera/data/gri30_highT.yaml", "gri30");

	Cantera::compositionMap initialComposition;
	initialComposition.emplace("O2",0.209);
	initialComposition.emplace("N2",1-0.209);

	freeStream.setTPX(300,101325,initialComposition);

	hypro::thermoState::warning_ = false;

	hypro::Node& cfdInlet(*new hypro::Node(freeStream));
	cfdInlet.Name_ = "CFD Inlet";
	hypro::Node& exit(*new hypro::Node(freeStream));
	exit.Name_ = "Nozz exit";

	std::unique_ptr<hypro::systemModule> sysPtr(new hypro::systemModule("SERJ",&cfdInlet,&exit));
	hypro::systemModule& system(*sysPtr);
	//system.verbosity_ = 0;

	double foot2 = std::pow(foot,2.0);
	//Nodes                            N1_                             0                             1
	const int numNodes = 9;
	double A[numNodes] =   {pi*(std::pow(1.14,2.0) - std::pow(0.01,2.0)),
			pi*((std::pow(1.14,2.0) - std::pow(0.29,2.0))-(std::pow(0.89,2.0) - std::pow(0.79,2.0))),
			pi*(std::pow(1.14,2.0) - std::pow(0.29,2.0)),pi*(std::pow(1.14,2.0) - std::pow(0.29,2.0)),
			0, pi*(std::pow(0.89,2.0) - std::pow(0.79,2.0)),
			pi*(std::pow(1.58,2.0) - std::pow(0.479,2.0)), 0.0, 0.0};
    //                  2         3             4                     5                6             N2_
	/*name = {'CFD Inlet','PinchPoint','Primary','Mixer End','Rocket Throat','Rocket Outlet',...
    'End Diffuser','Nozzle throat','Nozzle Exit'};*/

	cfdInlet.A(A[0]);
	system.initNodes(numNodes-2,cfdInlet);
	for(int i=0; i<numNodes-2; i++){
	    system.nodes_[i]->A(A[i+1]);
	}
	exit.A(A[numNodes-1]);

	double mfr, h;
	switch(selCase){
	case CASE1:
		//case 1 - h=0 M=0.75
		std::cout << "Selected Case: CASE1 h=0 M=0.75" << std::endl;
		// Set Inlet
		//cfdInlet.X(X);
		//cfdInlet.p_ = 26.12*psi;
		h = 136.12*btu/lbm;
        //cout << "h: " << h << endl;
		//cfdInlet.Th(h,600);
		cfdInlet.gasmodel_->setState_HP(h,26.12*psi);
		mfr = 1639*lbm;
		//cfdInlet.U_ = mfr/(cfdInlet.rho()*cfdInlet.A());
		cfdInlet.setU(mfr/(cfdInlet.rho()*cfdInlet.A()));

		// Set Nozzle areas
		system.nodes_[6]->A(66.44*foot2);
		exit.A(71.99*foot2);
		break;

	case CASE1tuned:
		//case 1 Tuned - h=0 M=0.75
		std::cout << "Selected Case: CASE1 tuned h=0 M=0.75" << std::endl;
		// Set Inlet
		//cfdInlet.X(X);
		//cfdInlet.p_ = 221273; //2.2315e+05;
		h = 136.2*btu/lbm;
		//cfdInlet.T_ =  610.251; //610.4113;
		std::cout << "CFD Inlet" << endl;
		cfdInlet.setTP(610.251,221273);
		mfr = 1639*lbm;
		//cfdInlet.U_ = mfr/(cfdInlet.rho()*cfdInlet.A());
		cfdInlet.setU(mfr/(cfdInlet.rho()*cfdInlet.A()));
		// Set Nozzle areas
		system.nodes_[6]->A(66.44*foot2);
		exit.A(71.99*foot2);
		break;

	case CASE2:
		//case 2
		std::cout << "Selected Case: CASE2 h=40000ft M=3.0" << std::endl;
		//set inlet
		//cfdInlet.X(X);
		//cfdInlet.p_ = 98.86*psi;
		h = 258.0*btu/lbm;
		//cfdInlet.Th(h,600);
		cfdInlet.gasmodel_->setState_HP(h,98.86*psi);
		mfr = 4552*lbm;
		//cfdInlet.U_ = mfr/(cfdInlet.rho()*cfdInlet.A());
		cfdInlet.setU(mfr/(cfdInlet.rho()*cfdInlet.A()));

		// Set Nozzle areas
		system.nodes_[6]->A(51.69*foot2);
		exit.A(106.93*foot2);
		break;
	}
/*
	//Print turbulence properties k and omega
	double k = 1.5*std::pow(cfdInlet.U_*0.1,2.0); //turbulent kinetic energy with 10% of turbulence intensity
	double omega = std::sqrt(k)/1; //turbulent specific dissipation rate with length scale=1m

	std::cout << "Inlet conditions:" << std::endl
			<< "\tgamma = " << cfdInlet.gamma() << std::endl
			<< "\tT0 = " << cfdInlet.gas().TH(cfdInlet.H(),600.0) << std::endl
			<< "\tflowRate = " << cfdInlet.mfr()*2/360 << std::endl
			<< "\tk = " << k << std::endl
			<< "\tomega = " << omega << std::endl;
	std::cout << "Parameters:" << std::endl << "\teta = " << eta << std::endl;
*/
/*
#ifdef MATLAB
	mexPrintf("Inlet conditions:\n\tflowRate = %g\n\tk = %g\n\tomega = %g\n",
	    cfdInlet.mfr()*2/360,k,omega);
#else
	std::printf("Inlet conditions:\n\tflowRate = %g\n\tk = %g\n\tomega = %g\n",
		    cfdInlet.mfr()*2/360,k,omega);
#endif
*/

	// Pre-primary
	hypro::systemModule::ModulePtr pinch(new hypro::isentropicDuct("pinch",&system.N1(),system.nodes_[0].get()));
	pinch->N2().Name_ = "Pinch Point";
	system.add(pinch);
	hypro::systemModule::ModulePtr prePrim(new hypro::isentropicDuct("pre-primary",system.nodes_[0].get(),system.nodes_[1].get()));
	prePrim->N2().Name_ = "Pre-Mixer";
	system.add(prePrim);

	// Primary - ejector
	hypro::systemModule::ModulePtr primaryNozzle(new hypro::isentropicDuct("primary",system.nodes_[3].get(),system.nodes_[4].get()));
	primaryNozzle->N1().Name_ = "Rocket-Throat";
	primaryNozzle->N2().Name_ = "Rocket-Out";
	mfr = 392*lbm;
	double p0 = 1500*psi;
	double T0 = 3675; //adiabatic flame calculation
	//primaryNozzle->N1().T_ = 1000; // Only to asses gamma
	//primaryNozzle->N1().X(cfdInlet.X());
	primaryNozzle->N1().setTPX(1000,101325,cfdInlet.X());
	double pstr_p0,Tstr_T0;
	hypro::isentropicDuct::isentropic(1.0,primaryNozzle->N1().gamma(),pstr_p0,Tstr_T0);
	//primaryNozzle->N1().p_ = p0*pstr_p0;
	//primaryNozzle->N1().T_ = T0*Tstr_T0;
	//cout << T0*Tstr_T0 << " " << p0*pstr_p0 << endl;
	primaryNozzle->N1().setTP(T0*Tstr_T0,p0*pstr_p0);
	//primaryNozzle->N1().U_ = primaryNozzle->N1().a();
	primaryNozzle->N1().setU(primaryNozzle->N1().geta());
	primaryNozzle->N1().A(mfr/(primaryNozzle->N1().rho()*primaryNozzle->N1().getU()));
	system.add(primaryNozzle);

	hypro::systemModule::ModulePtr primary(new hypro::Mixer("Mixer",system.nodes_[1].get(),
			system.nodes_[2].get(),system.nodes_[4].get()));
	primary->N2().Name_ = "Mixer-out";
	hypro::Mixer* primaryC = (hypro::Mixer*)primary.get();
	primaryC->teta_ = 0;
	primaryC->eff_ = eta;
	system.add(primary);

	// Diffuser
	hypro::systemModule::ModulePtr diffuser(new hypro::isentropicDuct("Diffuser",system.nodes_[2].get(),
			system.nodes_[5].get()));
	diffuser->N2().Name_ = "Diff-Out";
	system.add(diffuser);

	// Nozzle
	hypro::systemModule::ModulePtr nozzleConv(new hypro::chokedConv("Nozz-conv",system.nodes_[5].get(),
			system.nodes_[6].get(),hypro::chokedConv::CALCULATED));
	system.add(nozzleConv);
	system.nodes_[6]->Name_ = "Nozz Throat";
	hypro::systemModule::ModulePtr nozzle(new hypro::isentropicDuct("Nozzle",system.nodes_[5].get(),
			&exit));
	hypro::isentropicDuct* nozzleC = (hypro::isentropicDuct*)nozzle.get();
	nozzleC->choked_ = true;
	system.add(nozzle);

    system.verbosity_ = 0;
	system.calculate();

	return sysPtr;
}

void deleteAll(hypro::systemModule& sys){
	//const mixture< specieThermo< janafThermo<perfectGas> > >* gas = &sys.N1().specieGas();
	delete(&sys.N1());
	delete(&sys.N2());
	//delete(gas);
}
/*
#ifdef MATLAB

void mexFunction( int nlhs, mxArray *plhs[],
		int nrhs, const mxArray *prhs[] )
{
	size_t mrows,ncols;

	// Check for proper number of arguments.
	if(nrhs<1) {
		mexErrMsgIdAndTxt( "HYPRO:invalidNumInputs",
				"Not enough input arguments");
	}
	if(nlhs>3) {
		mexErrMsgIdAndTxt( "HYPRO:invalidNumOutputs",
				"Too many output arguments");
	}

	// The input must be a noncomplex scalar double.
        for(int i=0;i<nrhs;i++){
		mrows = mxGetM(prhs[i]);
		ncols = mxGetN(prhs[i]);
		if( !mxIsDouble(prhs[i]) || mxIsComplex(prhs[i]) ||
				!(mrows==1 && ncols==1) ) {
			mexErrMsgIdAndTxt( "HYPRO:inputNotRealScalarDouble",
					"Input must be a noncomplex scalar double.");
		}
	}

	//create matrix for the return argument.
	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);

	// Assign pointers to each input and output.
	int selCase = static_cast<int>(*mxGetPr(prhs[0]));
	double eta;
	if(nrhs>1){
		eta = *mxGetPr(prhs[1]);
	}else{
		eta = 1.0;
	}

	double& Thrust = *mxGetPr(plhs[0]);
	Foam::List<double> Mfr;
	Collection::OutType out;
	std::unique_ptr<systemModule> sys;

	try{
		feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
		sys = Run((caseTag)selCase,eta);
		Thrust = sys->thrust();
		Mfr = sys->propellantMfr();
		out = sys->out();
		fedisableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
	}catch (std::exception& e){
		fedisableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
		mexErrMsgIdAndTxt( "HYPRO:error", e.what());
	}
    
	if(nlhs>=2){
		plhs[1] = mxCreateDoubleMatrix((mwSize)mrows, Mfr.size(), mxREAL);
		double* point = mxGetPr(plhs[1]);
		for (mwSize i = 0; i < Mfr.size(); i++) {
			point[i] = Mfr[i];
		}
	}

	if(nlhs>=3){
		plhs[2] = mxCreateCellMatrix(out.size(), 2);
		for (mwIndex i = 0; i < (mwIndex)out.size(); i++){
			mwIndex I[2] = {i,0};
			mwIndex j = mxCalcSingleSubscript(plhs[2], 2, I);
			mxSetCell(plhs[2], j, mxCreateString(out[i].first.c_str()));

			I[1] = 1;
			j = mxCalcSingleSubscript(plhs[2], 2, I);
			mxArray* dat = mxCreateDoubleMatrix(1, out[i].second.size(), mxREAL);
			mxSetCell(plhs[2], j, dat);
			double* point = mxGetPr(dat);
			for (size_t k = 0; k < out[i].second.size(); k++) {
				point[k] = out[i].second[k];
			}
		}
	}

	deleteAll(*sys.get());
}

#else
*/

/*
int main(int argc, char *argv[]) {
	// Pharse inputs
	caseTag selCase = CASE1tuned;
	double eta = 1.0;
	for(int i=1; i<argc; i+=2){
		if(strcmp(argv[i],"-case")==0){
			sscanf(argv[i+1],"%d",(int*)&selCase);
		}else if(strcmp(argv[i],"-eta")==0){
			sscanf(argv[i+1],"%f",(float*)&eta);
		}else{
			std::cerr << "Error: wrong input." << std::endl;
			exit(1);
		}
	}

	std::unique_ptr<systemModule> sys = Run(selCase,eta);

	Foam::List<double> Mfr = sys->propellantMfr();
	Collection::OutType out(sys->out());

	sys->printOut(std::cout);

	deleteAll(*sys.get());
	return 1;
}

#endif
*/
