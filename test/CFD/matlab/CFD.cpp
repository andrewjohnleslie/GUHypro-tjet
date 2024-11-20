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

#include "CFD.h"
#include "mex.h"
#include <fenv.h>

void mexFunction( int nlhs, mxArray *plhs[],
		int nrhs, const mxArray *prhs[] )
{
	size_t mrows,ncols;

	/* Check for proper number of arguments. */
	if(nrhs<1) {
		mexErrMsgIdAndTxt( "HYPRO:invalidNumInputs",
				"Not enough input arguments");
	}
	if(nlhs>3) {
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
	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);

	/* Assign pointers to each input and output. */
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
		sys = CFD((caseTag)selCase,eta);
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
