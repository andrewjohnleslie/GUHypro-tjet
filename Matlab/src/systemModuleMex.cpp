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

#include <core/systemModule.h>
#include <fenv.h>
#include "mex.h"
#include <fstream>

#include <iostream>
#include <fstream>
#include <time.h>

#include "core/Node.h"
//#include "Collection.cpp"

#include <handleToPtr.h>

mxArray* constructor(std::string fileName, const mixture<specieThermo<janafThermo<perfectGas>>>& specieGas){
	std::ifstream fil(fileName);
	if(fil.fail()){
		std::string errStr("Error: cannot open file " + fileName);
		mexErrMsgIdAndTxt("HYPRO:invalidIO", errStr.c_str());
	}
	std::map<unsigned, std::shared_ptr<Node> > nodeMap;
	std::map<unsigned, std::shared_ptr<propModule> > moduleMap;
	systemModule::NodePtr N1 = Node::unserialize(fil,specieGas);
	systemModule::NodePtr N2 = Node::unserialize(fil,specieGas);
	nodeMap[N1->ID_] = N1;
	nodeMap[N2->ID_] = N2;

	Node* N1a = new Node(*N1.get());
	N1a->ID_ = N1->ID_;

	Node* N2a = new Node(*N2.get());
	N2a->ID_ = N2->ID_;


	systemModule::ModulePtr sys = propModule::unserialize(fil,nodeMap,moduleMap);
	fil.close();

	systemModule* sysa = new systemModule(*(systemModule*)sys.get());
	sysa->ID_ = sys->ID_;
	sysa->N1(*N1a);
	sysa->N2(*N2a);

	mxArray* handle  = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
	*reinterpret_cast<propModule**>(mxGetPr(handle)) = (propModule*)sysa;

	return handle;
}

void mexFunction( int nlhs, mxArray *plhs[],
		int nrhs, const mxArray *prhs[] )
{
	if (mxIsChar(prhs[0]) && (mxGetM(prhs[0]) == 1 ) ){
		size_t buflen = mxGetN(prhs[0])*sizeof(mxChar)+1;
		char* buf = (char*)mxMalloc(buflen);

		int status = mxGetString(prhs[0], buf, (mwSize)buflen);

		if(!strcmp(buf,"constructor")){
			if(nrhs!=3){
				mexErrMsgIdAndTxt( "HYPRO:invalidIO",
						"3 input required.");
			}
			if(nlhs>1) {
				mexErrMsgIdAndTxt( "HYPRO:invalidIO",
						"Too many output arguments");
			}

			if(!mxIsChar(prhs[1]) || (mxGetM(prhs[1]) != 1 )){
				mexErrMsgIdAndTxt( "HYPRO:invalidIO",
										"Second input must be a string.");
			}
			size_t strlen = mxGetN(prhs[1])*sizeof(mxChar)+1;
			char* fileName = (char*)mxMalloc(strlen);
			mxGetString(prhs[1], fileName, (mwSize)strlen);
			std::string fileNameStr(fileName);

			mixture<specieThermo<janafThermo<perfectGas>>>* specieGas =
					handleToPtr< mixture<specieThermo<janafThermo<perfectGas>>> >(prhs[2]);

			try{
				feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
				plhs[0] = constructor(fileNameStr,*specieGas);
				fedisableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
			}catch(std::exception& e){
				fedisableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
				mexErrMsgIdAndTxt( "HYPRO:error", e.what());
			}

			return;

		}else if(!strcmp(buf,"destructor")){
				mexErrMsgIdAndTxt( "HYPRO:invalidIO",
						"Destructor not implemented.");
			return;

		}else if(!strcmp(buf,"thrust")){
			if(nrhs!=2){
				mexErrMsgIdAndTxt( "HYPRO:invalidNumOutputs",
						"2 input required.");
			}
			if(nlhs>1) {
				mexErrMsgIdAndTxt( "HYPRO:invalidNumOutputs",
						"Too many output arguments");
			}
			plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
			double& out = *mxGetPr(plhs[0]);
			try{
				feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
				systemModule* obj = (systemModule*)handleToPtr<propModule>(prhs[1]);
				out = obj->thrust();
				fedisableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
			}catch(std::exception& e){
				fedisableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
				mexErrMsgIdAndTxt( "HYPRO:error", e.what());
			}
			return;

		}else if(!strcmp(buf,"propellantMfr")){
			if(nrhs!=2){
				mexErrMsgIdAndTxt( "HYPRO:invalidNumOutputs",
						"2 input required.");
			}
			if(nlhs>1) {
				mexErrMsgIdAndTxt( "HYPRO:invalidNumOutputs",
						"Too many output arguments");
			}

			Foam::List<double> out;
			try{
				systemModule* obj = (systemModule*)handleToPtr<propModule>(prhs[1]);
				out = obj->propellantMfr();
			}catch(std::exception& e){
				mexErrMsgIdAndTxt( "HYPRO:error", e.what());
			}

			plhs[0] = mxCreateDoubleMatrix(1, out.size(), mxREAL);
			double* point = mxGetPr(plhs[0]);
			for (mwSize i = 0; i < out.size(); i++) {
				point[i] = out[i];
			}
			return;

		}else if(!strcmp(buf,"nodes")){
			if(nrhs!=2){
				mexErrMsgIdAndTxt( "HYPRO:invalidNumOutputs",
						"2 input required.");
			}
			if(nlhs>1) {
				mexErrMsgIdAndTxt( "HYPRO:invalidNumOutputs",
						"Too many output arguments");
			}

			systemModule* obj = (systemModule*)handleToPtr<propModule>(prhs[1]);

			plhs[0]  = mxCreateNumericMatrix(obj->nodes_.size(), 1, mxUINT64_CLASS, mxREAL);
			for (std::size_t i = 0; i < obj->nodes_.size(); i++) {
				*reinterpret_cast<thermoState**>(&(mxGetPr(plhs[0])[i])) = (thermoState*)obj->nodes_[i].get();
			}
			return;
		}else if(!strcmp(buf,"modules")){
					if(nrhs!=2){
						mexErrMsgIdAndTxt( "HYPRO:invalidNumOutputs",
								"2 input required.");
					}
					if(nlhs>1) {
						mexErrMsgIdAndTxt( "HYPRO:invalidNumOutputs",
								"Too many output arguments");
					}

					systemModule* obj = (systemModule*)handleToPtr<propModule>(prhs[1]);

					plhs[0]  = mxCreateNumericMatrix(obj->modules_.size(), 1, mxUINT64_CLASS, mxREAL);
					//for (std::size_t i = 0; i < obj->modules_.size(); i++) {
					//	*reinterpret_cast<propModule**>(&(mxGetPr(plhs[0])[i])) = (propModule*)obj->modules_[i].get();
					//}
					int i = 0;
					for(std::list<systemModule::ModulePtr>::iterator it=obj->modules_.begin(); it!=obj->modules_.end(); it++){
						*reinterpret_cast<propModule**>(&(mxGetPr(plhs[0])[i])) = it->get();
						i++;
					}
			return;
		}else{
			std::ostringstream er;
			er << "Flag " << buf << " does not exist.";
			mexErrMsgIdAndTxt( "HYPRO:invalidNumInputs", er.str().c_str());
		}

		mxFree(buf);
	}else{
		mexErrMsgIdAndTxt( "HYPRO:invalidInputs", "The first input must be a string.");
	}
}
