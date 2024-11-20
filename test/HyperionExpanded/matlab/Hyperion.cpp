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

#include "Hyperion.h"

#include "mex.h"
#include <fenv.h>


void mexFunction( int nlhs, mxArray *plhs[],
		int nrhs, const mxArray *prhs[] )
{
	size_t mrows,ncols;
    std::cout << nlhs << " " << nrhs << std::endl;
    std::cout << "Here1" << std::endl;
	/* Check for proper number of arguments. */
	if(nrhs!=6) {
		mexErrMsgIdAndTxt( "HYPRO:invalidNumInputs",
				"Number of inputs required: 6");
	}
	if(nlhs>6) {
		mexErrMsgIdAndTxt( "HYPRO:invalidNumOutputs",
				"Too many output arguments");
	}
    std::cout << "Here2" << std::endl;

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
    std::cout << "Here3" << std::endl;
	/* Create matrix for the return argument. */
	//plhs[0] = mxCreateDoubleMatrix((mwSize)mrows, (mwSize)ncols, mxREAL);
	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	plhs[3] = mxCreateDoubleMatrix(1, 1, mxREAL);
	plhs[4] = mxCreateDoubleMatrix(1, 1, mxREAL);

	/* Assign pointers to each input and output. */
	std::cout << "Here3.1" << std::endl;
	const double p = *mxGetPr(prhs[0]);
	std::cout << "Here3.2, p = " << p << std::endl;
	const double T = *mxGetPr(prhs[1]);
	std::cout << "Here3.3, T = " << T << std::endl;
	const double Mach = *mxGetPr(prhs[2]);
	std::cout << "Here3.4, Mach = " << Mach << std::endl;
    const double mode = *mxGetPr(prhs[3]);
	std::cout << "Here3.5, Mode = " << mode << std::endl;
    const double throttle = *mxGetPr(prhs[4]);
	std::cout << "Here3.6, throttle = " << throttle << std::endl;
    const double verbosity = *mxGetPr(prhs[5]);
    std::cout << "Here4, verbosity = " << verbosity << std::endl;
//	int mode = -1;
//	if(nrhs==4){
//		double mode1 = *mxGetPr(prhs[3]);
//		mode1 += 0.5;
//		mode = (int)mode1;
//	}
	double& Thrust = *mxGetPr(plhs[0]);
//	Foam::List<double> Mfr;
    std::vector<double> Mfr;
    std::cout << "Here4.1" << std::endl;
	Collection::OutType out;
    std::cout << "Here4.2" << std::endl;
    double& dynamic_pressure = *mxGetPr(plhs[3]);
    std::cout << "Here4.3" << std::endl;
    double& total_flow = *mxGetPr(plhs[4]);
    std::cout << "Here4.4" << std::endl;
    std::vector<double> emissions;
    std::cout << "Here5" << std::endl;

	systemModule* sys;
	try{
    //    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
        systemModule *sys = Hyperion(p, T, Mach, throttle, mode, verbosity, Thrust, Mfr, out, dynamic_pressure, emissions, total_flow);
     //   fedisableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

    }catch (std::exception& e){
		fedisableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
		mexErrMsgIdAndTxt( "HYPRO:error", e.what());
	}

	std::cout << "Here 6" << std::endl;

    // MFR
    plhs[1] = mxCreateDoubleMatrix((mwSize)mrows, Mfr.size(), mxREAL);
	double* point = mxGetPr(plhs[1]);
	for (mwSize i = 0; i < Mfr.size(); i++) {
		point[i] = Mfr[i];
    }


    // Collection out
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

    plhs[5] = mxCreateDoubleMatrix((mwSize)mrows, emissions.size(), mxREAL);
    point = mxGetPr(plhs[5]);
    for( mwSize i = 0; i < emissions.size(); i++){
        point[i] = Mfr[i];
    }


//    if(nlhs>=7){
//    	plhs[6]  = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
//    	*reinterpret_cast<propModule**>(mxGetPr(plhs[3])) = (propModule*)sys;
//    }else{
//    	systemModule::deleteAll(sys);
//    }
}
