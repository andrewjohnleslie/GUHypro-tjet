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

#include "SCRAMSPACE.h"
#include "mex.h"
#include <fenv.h>

void mexFunction( int nlhs, mxArray *plhs[],
		int nrhs, const mxArray *prhs[] )
{
	/* Check for proper number of arguments. */
	if(nlhs>1) {
		mexErrMsgIdAndTxt( "HYPRO:invalidNumOutputs",
				"Too many output arguments");
	}
	if(nrhs<1){
		mexErrMsgIdAndTxt( "HYPRO:invalidIO",
						"Not enough input arguments");
	}

	/* Assign pointers to each input and output. */
	const double p = 4100;
	const double T = 370;
	const double U = 2830;
	const double combEff = *mxGetPr(prhs[0]);

	double Thrust;
	Foam::List<double> Mfr;
	Collection::OutType out;

	try{
		feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
		SCRAMSPACE(p,T,U,combEff,0,Thrust,Mfr,out);
		fedisableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
	}catch (std::exception& e){
		fedisableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
		mexErrMsgIdAndTxt( "HYPRO:error", e.what());
	}

    if(nlhs>=1){
    	plhs[0] = mxCreateCellMatrix(out.size(), 2);
    	for (mwIndex i = 0; i < (mwIndex)out.size(); i++){
    		mwIndex I[2] = {i,0};
    		mwIndex j = mxCalcSingleSubscript(plhs[0], 2, I);
    		mxSetCell(plhs[0], j, mxCreateString(out[i].first.c_str()));

    		I[1] = 1;
    		j = mxCalcSingleSubscript(plhs[0], 2, I);
    		mxArray* dat = mxCreateDoubleMatrix(1, out[i].second.size(), mxREAL);
    		mxSetCell(plhs[0], j, dat);
    		double* point = mxGetPr(dat);
    		for (size_t k = 0; k < out[i].second.size(); k++) {
    			point[k] = out[i].second[k];
    		}
    	}
    }
}
