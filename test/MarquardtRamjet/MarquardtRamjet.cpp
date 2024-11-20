/*!
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

#define WM_DP
#define NoRepository

#define DO_QUOTE(X)       #X
#define QUOTE(X)          DO_QUOTE(X)

#include <iostream>
#include <fstream>
#include <time.h>

#include "core/Node.h"

#include <IOdictionary.H>
#include <IOobject.H>
#include <Time.H>
#include <argList.H>
#include <specieThermo.H>
#include <janafThermo.H>
#include <perfectGas.H>
#include <perfectGas.C>
using namespace Foam;

#include <mixture.h>
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
#include "solvers/EffIsenDuct.h"
#include "solvers/Friction.h"
#include "ISAModel.h"

#ifdef MATLAB
#include "mex.h"
#include <fenv.h>
#endif

systemModule* Hyperion(const double& p, const double& T, const double& Mach, const int mode,
		double& Thrust, Foam::List<double>& Mfr, Collection::OutType& out){

	double foot =  0.3048; //conversion feet to meters
	double lbm = 0.45359; //conversion libre to kg

	// Thermal model

	//Time runTime("/home/trb12187/Matlab/HyPro","thermo");
	Time runTime(QUOTE(HyProRescource),"");

	IOdictionary a
	(
			IOobject
			(
					"thermophysicalProperties",
					"",
					runTime,
					IOobject::MUST_READ,
					IOobject::NO_WRITE,
					false
			)
	);
	wordList species(4,"H2");
	species[1] = "O2";
	species[2] = "N2";
	species[3] = "H2O";

	mixture< specieThermo< janafThermo<perfectGas> > > specieGas(a,species);

	Node freeStream(specieGas);
	freeStream.p_ = p;
	freeStream.T_ = T;
	List<double> X(4,0.0);
	X[1] = 0.209;   //Air
	X[2] = 1-0.209;
	freeStream.X(X);
	freeStream.U_ = Mach*freeStream.a();
	freeStream.Name_ = "Free Stream";

	thermoState::warning_ = false;

	Node freeStreamMod(specieGas);
	freeStreamMod.Name_ = "Pre Intake";
	Node exit(specieGas);
	exit.Name_ = "Nozzle Exit";

	Wedge forebody("forebody",&freeStream,&freeStreamMod);
	forebody.delta_ = 8.0;
	forebody.calculate();

	MultiModeMachP0 multiSystem("system",&freeStreamMod,&exit);
	multiSystem.verbosity_ = 1;
	systemModule& system(multiSystem.add("Ejector Mode"));
	system.verbosity_ = 0;

	//Nodes                 N1_    I0    I1       0     1       2     3      4            5       6     7     N2_
	const int numNodes = 12;
	double A[numNodes] =   {70.0, 70.0, 57.0, 	27.02, 43.49, 43.49, 0.0, (43.49-27.02), 76.24, 76.24, 76.24, 107.0};
	/*name = {'Pre Intake','Intake','Throat','Pinch Point','Primary','Mixer End',...
	    'Dummy','Rocket Outlet','End Diffuser','Injection','End Chamber','Nozzle'};*/
	double foot2 = std::pow(foot,2.0);
	freeStreamMod.A(A[0]*foot2);
	system.initNodes(8,freeStreamMod);
	for(int i=0; i<8; i++){
	    system.nodes_[i]->A(A[i+3]*foot2);
	}
	exit.A(A[numNodes-1]*foot2);
	exit.Amax_ = A[numNodes-1]*foot2;

	//Inlet
	systemModule::ModulePtr inlet(new AdaptedThroatInlet("Inlet",&system.N1(),system.nodes_[0].get()));
	inlet->N2().Name_ = "Pinch Point";
	AdaptedThroatInlet* inletC = (AdaptedThroatInlet*)inlet.get();
	inletC->nodes_[0]->A(A[1]*foot2);
	inletC->nodes_[0]->Amax_ = inletC->nodes_[0]->A();
	inletC->nodes_[0]->Name_ = "Intake";
	inletC->nodes_[1]->A(A[2]*foot2);
	inletC->nodes_[1]->Amax_ = inletC->nodes_[1]->A();
	inletC->nodes_[1]->Amin_ = 0;
	inletC->nodes_[1]->Name_ = "Throat";
	inlet->verbosity_ =0;
	inlet->chokingFeedback_ = inlet;
	system.add(inlet);

	// Primary - ejector
	systemModule::ModulePtr postPinch(new isentropicDuct(
			"Post pinch point",system.nodes_[0].get(),system.nodes_[1].get()));
	system.add(postPinch);

	systemModule::ModulePtr rocket(new Rocket("Rocket",system.nodes_[3].get(),system.nodes_[4].get()));
	rocket->N2().Name_ = "Rocket exit";
	Rocket* rocketC = (Rocket*)rocket.get();
	List<double> Xr(4,0.0);
	Xr[3] = 1;
	rocketC->N2().X(Xr);
	rocketC->Me_ = 6.0;
	rocketC->mfr_ = 216*lbm;
	rocketC->Isp_ = 462;
	system.add(rocket);

	systemModule::ModulePtr primary(new Mixer("primary",
			system.nodes_[1].get(),system.nodes_[2].get(),system.nodes_[4].get()));
	primary->N1().Name_ = "Primary In";
	primary->N2().Name_ = "Primary Out";
	Mixer* primaryC = (Mixer*)primary.get();
	primaryC->verbosity_ = 0;
	primaryC->chokingFeedback_ = inlet;
	primaryC->eff_ = 0.85;
	system.add(primary);

	// Secondary - Post combustor
	systemModule::ModulePtr diffuser(new EffIsenDuct("diffuser",system.nodes_[2].get(),system.nodes_[5].get()));
	diffuser->N2().Name_ = "Diffuser End";
	EffIsenDuct* diffuserC = (EffIsenDuct*)diffuser.get();
	diffuserC->etap_ = 0.7;
	diffuserC->etaT_ = 1.0;
	system.add(diffuser);

	systemModule::ModulePtr secondary(new InjectionPhi("secondary",system.nodes_[5].get(),system.nodes_[6].get()));
	secondary->N2().Name_ = "Injection";
	InjectionPhi* secondaryC = (InjectionPhi*)secondary.get();
	secondaryC->PhiMax_ = 1;
	secondaryC->Phi_ = 1;
	List<double> Xf(4,0.0);
	Xf[0] = 1;
	secondaryC->Xf_ = Xf;
	secondaryC->Tf_ = 20;
	secondary->chokingFeedback_ = secondary;
	system.add(secondary);

	systemModule::ModulePtr combustor(new Combustion<Friction>("combustor",system.nodes_[6].get(),system.nodes_[7].get()));
	combustor->N2().Name_ = "Chamber End";
	combustor->chokingFeedback_ = secondary;
	Combustion<Friction>* combustorC = (Combustion<Friction>*)combustor.get();
	combustorC->eta_ = 0.9;
	combustorC->Cf_ = 0.01;
	combustorC->D_ = 1.631;
	combustorC->L_ = 1.0;
	system.add(combustor);

	// Nozzle
	systemModule::ModulePtr nozzle(new AdaptedNozzle<EffIsenDuct>("nozzle",system.nodes_[7].get(),&system.N2()));
	AdaptedNozzle<EffIsenDuct>* nozzleC = (AdaptedNozzle<EffIsenDuct>*)nozzle.get();
	nozzleC->choked_ = true;
	/* in case of inlet spilling (especially in supersonic) better to adapt the
		nozzle to the inlet conditions.*/
	nozzleC->freeStream_ = &inlet->N1();
	nozzleC->etap_ = 1.0;
	system.add(nozzle);

	// Ramjet mode
	systemModule& ramjet(multiSystem.add("Ramjet Mode",0));
	//const systemModule::ModulePtr inlet1(ramjet.exchange<AdaptedInlet>(inlet,"Inlet"));
	//AdaptedInlet* inlet1C = (AdaptedInlet*)inlet1.get();
	//*(inlet1C->nodes_[0]) = *(inletC->nodes_[0]);
	//inlet1C->nodes_[1]->A(inlet1C->nodes_[1]->Amax_);
	//inlet1C->chokingFeedback_ = inlet1;

	const systemModule::ModulePtr primary1(ramjet.exchange<NeutralLink>(primary,"primary"));
	ramjet.remove(rocket);

	//Pure Rocket mode
	systemModule& rocketMode(multiSystem.add("Rocket Mode",0));
	const systemModule::ModulePtr primary2(rocketMode.exchange<isentropicDuct>(primary,"primary"));
	primary2->N1(rocket->N2());
	const systemModule::ModulePtr secondary2(rocketMode.exchange<NeutralLink>(secondary,"secondary"));
	secondary2->N2(combustor->N2());
	rocketMode.remove(inlet);
	rocketMode.remove(postPinch);
	rocketMode.remove(combustor);
	Node rocketDummyN1(rocketMode.N1());
	rocketMode.N1(rocketDummyN1);
	rocketDummyN1.A(0.0);

	multiSystem.rangesM_[0][0] = 0.0;
	multiSystem.rangesM_[0][1] = 3.0;
	multiSystem.rangesM_[1][0] = 3.0;
	multiSystem.rangesM_[1][1] = 7.0;
	multiSystem.rangesP0_[0][0] = 1e4;
	multiSystem.rangesP0_[1][0] = 5e5;


	multiSystem.unreduce();

	if(mode==-1){
		multiSystem.calculate();
		Thrust = multiSystem.thrust();
		Mfr = multiSystem.propellantMfr();
	}else{
		systemModule& selMode(multiSystem.modes_[mode]);
		selMode.verbosity_ = 1;
		selMode.calculate();
		Thrust = selMode.thrust();
		Mfr = selMode.propellantMfr();
	 	system.printOut(std::cout);

		std::array<double,10> outTmp;
		outTmp[0] = freeStream.A();
		outTmp[1] = freeStream.Amax_;
		outTmp[2] = freeStream.Amin_;
		outTmp[3] = freeStream.p_;
		outTmp[4] = freeStream.T_;
		outTmp[5] = freeStream.U_;
		outTmp[6] = freeStream.X()[0];
		outTmp[7] = freeStream.X()[1];
		outTmp[8] = freeStream.X()[2];
		outTmp[9] = freeStream.X()[3];
		out.push_back(std::make_pair(freeStream.Name_,outTmp));


		/*std::cout << "T0:" << std::endl;

		for(std::size_t i=0; i<system.nodes_.size(); i++){
		double H = system.nodes_[i]->H();
		double T0 = system.nodes_[i]->gas().TH(H,300);
		std::cout << system.nodes_[i]->Name_ << " " << T0 << std::endl;

			}

		double H = system.N2().H();
		double T0 = system.N2().gas().TH(H,300);
		std::cout << "Nozzle Exit " << T0 << std::endl;*/


		Collection::OutType out1(selMode.out());
		for(std::size_t i=0; i<out1.size(); i++){
			out.push_back(out1[i]);
		}
	}

	/*std::ofstream fil("Hyperion.hyp");
	system.N1().serialize(fil);
	fil << std::endl;
	system.N2().serialize(fil);
	fil << std::endl;
	system.serialize(fil);
	fil.close();*/

	/*std::ofstream fil("HyperionRamjet.hyp");
	ramjet.N1().serialize(fil);
	fil << std::endl;
	ramjet.N2().serialize(fil);
	fil << std::endl;
	ramjet.serialize(fil);
	fil.close();*/

	return new systemModule(system);
}

#ifdef MATLAB

void mexFunction( int nlhs, mxArray *plhs[],
		int nrhs, const mxArray *prhs[] )
{
	size_t mrows,ncols;

	/* Check for proper number of arguments. */
	if(nrhs<3) {
		mexErrMsgIdAndTxt( "HYPRO:invalidNumInputs",
				"Not enough input arguments");
	}
	if(nlhs>4) {
		mexErrMsgIdAndTxt( "HYPRO:invalidNumOutputs",
				"Too many output arguments");
	}

	/* The input must be a noncomplex scalar double.*/
	for(int i=0;i<nrhs;i++){
		mrows = mxGetM(prhs[i]);
		ncols = mxGetN(prhs[i]);
		if( !mxIsDouble(prhs[i]) || mxIsComplex(prhs[i]) ||
				!(mrows==1 && ncols==1) ) {
			mexErrMsgIdAndTxt( "HYPRO:inputNotRealScalarDouble",
					"Input must be a noncomplex scalar double.");
		}
	}

	/* Create matrix for the return argument. */
	plhs[0] = mxCreateDoubleMatrix((mwSize)mrows, (mwSize)ncols, mxREAL);

	/* Assign pointers to each input and output. */
	const double p = *mxGetPr(prhs[0]);
	const double T = *mxGetPr(prhs[1]);
	const double Mach = *mxGetPr(prhs[2]);
	int mode = -1;
	if(nrhs==4){
		mode = *(int*)mxGetPr(prhs[3]);
	}
	double& Thrust = *mxGetPr(plhs[0]);
	Foam::List<double> Mfr;
	Collection::OutType out;

	systemModule* sys;
	try{
		feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
		sys = Hyperion(p,T,Mach,mode,Thrust,Mfr,out);
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

    if(nlhs>=4){
    	plhs[3]  = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
    	*reinterpret_cast<propModule**>(mxGetPr(plhs[3])) = (propModule*)sys;
    }
}

#else

int main(int argc, char *argv[]) {

	int Machiterations = 9;
	int Altiterations = 1;

	double MachRange[9] = {2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0};
	//double AltitudeRange[4] = {9144, 15240, 21336, 33258};

	//double MachRange[1] = {3.5};
	double AltitudeRange[1] = {9144};

	ISAModel ISA;


 	for(int j=0; j<Altiterations; j++){
 		double Altitude_ = AltitudeRange[j];
 		ISA.calculate(Altitude_);

 		double p = ISA.Pressure(Altitude_);
	 	double T = ISA.Temperature(Altitude_);

	 	std::cout  << std::endl << std::endl << "ALTITUDE: " << Altitude_ << std::endl << "Pressure: " << p <<" Temperature: " << T << std::endl << "Local Speed of Sound: " << std::pow(1.4*287*T,0.5);



	 	for(int i=0; i<Machiterations; i++){

	 		double Mach = MachRange[i];

	 		std::cout  << std::endl << "NEW" << std::endl <<  "Mach Number: " << Mach << std::endl << std::endl;

	 		double Thrust;
	 		Foam::List<double> Mfr;
	 		Collection::OutType out;

	 		Hyperion(p,T,Mach,1,Thrust,Mfr,out);


	 		std::cout  << std::endl <<"Thrust = " << Thrust << ", Specific Impulse = " << Thrust/(9.81*Mfr[0]) <<  std::endl;
	 		//std::cout  << std::endl << Thrust;

		}

 	}

 	return 1;
}


#endif
