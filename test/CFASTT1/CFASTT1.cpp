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
#include <ClosedInlet.h>

#ifdef MATLAB
#include "mex.h"
#include <fenv.h>
#endif

MultiMode* Create(){
	thermoState::warning_ = false;

	double foot =  0.3048; //conversion feet to meters
	double lbm = 0.45359; //conversion libre to kg
	double inch = 0.0254;

	// Thermal model

	std::cout << QUOTE(HyProSrcPth) << std::endl;
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

	mixture< specieThermo< janafThermo<perfectGas> > >* specieGas =
			new mixture< specieThermo< janafThermo<perfectGas> > >(a,species);

	Node* freeStream = new Node(*specieGas);
	freeStream->p_ = 0.0;
	freeStream->T_ = 0.0;
	List<double> X(4,0.0);
	X[1] = 0.209;   //Air
	X[2] = 1-0.209;
	freeStream->X(X);
	freeStream->U_ = 0.0;

	Node* nozzExit = new Node(*specieGas);

	MultiMode* multiSystem = new MultiModeMachP0("system",freeStream,nozzExit);
	multiSystem->verbosity_ = 0;
	systemModule& system(multiSystem->add("Ejector Mode"));
	system.verbosity_ = 0;

	//double K = 0.3;
	double Aprim = 11.25 - 8.24;
	double Asec = 27;
	//Nodes                 N1_    I0    I1       0     1       2     3    4    5      6     7   N2_
	const int numNodes = 12;
	double A[numNodes] =   {27.0, 27.0, 0.25*27, 11.25, 11.25, 11.25, 0.0, 0, 22.5, 22.5, 22.5, 95.0};
	/*name = {'Pre Intake','Intake','Throat','Pinch Point','Primary','Mixer End',...
		    'Dummy','Rocket Outlet','End Diffuser','Injection','End Chamber','Nozzle'};*/
	for(int i=0; i<numNodes; i++){
		A[i] = A[i];
	}
	double foot2 = std::pow(foot,2.0);
	A[3] = A[3] - Aprim;
	A[7] = Aprim;

	freeStream->A(A[0]*foot2);
	system.initNodes(8,*freeStream);
	for(int i=0; i<8; i++){
		system.nodes_[i]->A(A[i+3]*foot2);
	}
	nozzExit->A(A[numNodes-1]*foot2);
	nozzExit->Amax_ = A[numNodes-1]*foot2;

	//Inlet
	systemModule::ModulePtr inlet(new AdaptedThroatInlet("Inlet",&system.N1(),system.nodes_[0].get()));
	AdaptedThroatInlet* inletC = (AdaptedThroatInlet*)inlet.get();
	inletC->nodes_[0]->A(A[1]*foot2);
	inletC->nodes_[0]->Amax_ = inletC->nodes_[0]->A();
	inletC->nodes_[1]->A(A[2]*foot2);
	inletC->nodes_[1]->Amax_ = inletC->nodes_[1]->A();
	inletC->chokingFeedback_ = inlet;
	inlet->verbosity_ =0;
	system.add(inlet);

	// Primary - ejector
	systemModule::ModulePtr postPinch(new isentropicDuct(
			"Post pinch point",system.nodes_[0].get(),system.nodes_[1].get()));
	system.add(postPinch);

	systemModule::ModulePtr rocket(new Rocket("Rocket",system.nodes_[3].get(),system.nodes_[4].get()));
	Rocket* rocketC = (Rocket*)rocket.get();
	List<double> Xr(4,0.0);
	Xr[3] = 1;
	rocketC->N2().X(Xr);
	rocketC->Me_ = 6.0;
	rocketC->mfrMax_ = 216*lbm;
	rocketC->Isp_ = 462;
	rocketC->unreduce();
	system.add(rocket);

	systemModule::ModulePtr primary(new Mixer("primary",
			system.nodes_[1].get(),system.nodes_[2].get(),system.nodes_[4].get()));
	Mixer* primaryC = (Mixer*)primary.get();
	primaryC->verbosity_ = 0;
	primaryC->chokingFeedback_ = inlet;
	primaryC->eff_ = 0.67;
	system.add(primary);

	// Secondary - Post combustor
	systemModule::ModulePtr diffuser(new isentropicDuct("diffuser",system.nodes_[2].get(),system.nodes_[5].get()));
	system.add(diffuser);

	systemModule::ModulePtr secondary(new InjectionPhi("secondary",system.nodes_[5].get(),system.nodes_[6].get()));
	InjectionPhi* secondaryC = (InjectionPhi*)secondary.get();
	secondaryC->PhiMax_ = 1;
	secondaryC->Phi_ = 1;
	List<double> Xf(4,0.0);
	Xf[0] = 1;
	secondaryC->Xf_ = Xf;
	secondary->chokingFeedback_ = secondary;
	system.add(secondary);

	systemModule::ModulePtr combustor(new Combustion<balanceMach>("combustor",system.nodes_[6].get(),system.nodes_[7].get()));
	combustor->chokingFeedback_ = secondary;
	system.add(combustor);

	// Nozzle
	systemModule::ModulePtr nozzle(new AdaptedNozzle<isentropicDuct>("nozzle",system.nodes_[7].get(),&system.N2()));
	AdaptedNozzle<isentropicDuct>* nozzleC = (AdaptedNozzle<isentropicDuct>*)nozzle.get();
	nozzleC->choked_ = true;
	//nozzle.p0Ratio(1);
	/*in case of inlet spilling (especially in supersonic) better to adapt the
			nozzle to the inlet conditions.*/
	nozzleC->freeStream_ = &inlet->N1();
	system.add(nozzle);

	// Ramjet mode
	systemModule& ramjet(multiSystem->add("Ramjet Mode",0));
	const systemModule::ModulePtr inlet1(ramjet.exchange<AdaptedInlet>(inlet,"Inlet"));
	AdaptedInlet* inlet1C = (AdaptedInlet*)inlet1.get();
	*(inlet1C->nodes_[0]) = *(inletC->nodes_[0]);
	inlet1C->nodes_[1]->Amax_ = inlet1C->nodes_[0]->Amax_;
	inlet1C->nodes_[1]->A(inlet1C->nodes_[1]->Amax_);
	inlet1C->chokingFeedback_ = inlet1;

	const systemModule::ModulePtr primary1(ramjet.exchange<NeutralLink>(primary,"primary"));
	ramjet.remove(rocket);

	//Pure Rocket mode
	systemModule& rocketMode(multiSystem->add("Rocket Mode",0));
	const systemModule::ModulePtr primary2(rocketMode.exchange<isentropicDuct>(primary,"primary"));
	primary2->N1(rocket->N2());
	const systemModule::ModulePtr secondary2(rocketMode.exchange<NeutralLink>(secondary,"secondary"));
	secondary2->N2(combustor->N2());
	const systemModule::ModulePtr inlet2(rocketMode.exchange<ClosedInlet>(inlet,"inlet"));
	rocketMode.remove(postPinch);
	rocketMode.remove(combustor);
	const systemModule::ModulePtr nozzle2(rocketMode.exchange<AdaptedNozzle<isentropicDuct> >(nozzle,"nozzle"));
	AdaptedNozzle<isentropicDuct>* nozzle2C = (AdaptedNozzle<isentropicDuct>*)nozzle2.get();
	nozzle2C->choked_ = false;
	nozzle2C->freeStream_ = &multiSystem->N1();

	//Off Mode
	systemModule& off(multiSystem->add("Off Mode"));
	systemModule::ModulePtr link(new NeutralLink("Off",&off.N1(),&off.N2()));
	off.add(link);

	MultiModeMachP0* multiSystemC = (MultiModeMachP0*)multiSystem;
	multiSystemC->rangesM_[0][0] = 0.0;
	multiSystemC->rangesM_[0][1] = 3.0;
	multiSystemC->rangesM_[1][0] = 2.5;
	multiSystemC->rangesM_[1][1] = 8.0;
	multiSystemC->rangesP0_[0][0] = 1e4;
	multiSystemC->rangesP0_[1][0] = 5e5;

	return multiSystem;
}

#ifdef MATLAB
const double& checkGet(const mxArray* input){
	/* The input must be a noncomplex scalar double.*/

	size_t mrows = mxGetM(input);
	size_t ncols = mxGetN(input);
	if( !mxIsDouble(input) || mxIsComplex(input) ||
			!(mrows==1 && ncols==1) ) {
		mexErrMsgIdAndTxt( "HYPRO:inputNotRealScalarDouble",
				"Input must be a non-complex scalar double.");
	}

	/* Assign pointers to each input and output. */
	const double& out = *mxGetPr(input);
	if(mxIsNaN(out)){
		mexErrMsgIdAndTxt( "HYPRO:inputNotRealScalarDouble",
				"Input must not be NaN.");
	}
	return out;
}

void mexFunction( int nlhs, mxArray *plhs[],
		int nrhs, const mxArray *prhs[] )
{
	mexLock();

	if(nrhs==0) {
		mexErrMsgIdAndTxt( "HYPRO:invalidNumInputs",
				"Not enough input parameters");
	}

	/* Check for proper input type */
	if (mxIsChar(prhs[0]) && (mxGetM(prhs[0]) == 1 ) )  {
		size_t buflen = mxGetN(prhs[0])*sizeof(mxChar)+1;
		char* buf = (char*)mxMalloc(buflen);

		int status = mxGetString(prhs[0], buf, (mwSize)buflen);

		if(!strcmp(buf,"Create")){
			if(nrhs!=1){
				mexErrMsgIdAndTxt( "HYPRO:invalidNumOutputs",
						"1 input required.");
			}
			if(nlhs>1){
				mexErrMsgIdAndTxt( "HYPRO:invalidNumOutputs",
						"Too many output arguments");
			}
			if(nlhs==1){
				plhs[0]  = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
				MultiMode* multiSystem = Create();
				*reinterpret_cast<MultiMode**>(mxGetPr(plhs[0])) = multiSystem;
			}
			return;

		}else{
			std::ostringstream er;
			er << "Flag " << buf << " does not exist.";
			mexErrMsgIdAndTxt( "HYPRO:invalidNumInputs", er.str().c_str());
		}

		mxFree(buf);
	}
}
#endif
