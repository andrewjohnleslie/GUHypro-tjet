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
#include "US76.h"

//This is an excerpt (modified) from PDAS, the original code can be found
//online at http://www.pdas.com/atmosdownload.html

US76::US76(double alt) {
    const double REARTH=6369.0;    // radius of the Earth (km)
    const double GMR = 34.163195;
    const int NTAB = 8;
    int i,j,k;

    static double htab[NTAB] = {0.0,  11.0, 20.0, 32.0, 47.0,
                                51.0, 71.0, 84.852 };
    static double ttab[NTAB] = { 288.15, 216.65, 216.65, 228.65, 270.65,
                                 270.65, 214.65, 186.946 };
    static double ptab[NTAB] = { 1.0, 2.2336110E-1, 5.4032950E-2, 8.5666784E-3,
                                 1.0945601E-3, 6.6063531E-4, 3.9046834E-5, 3.68501E-6 };
    static double gtab[NTAB] = { -6.5, 0, 1.0, 2.8, 0, -2.8, -2.0, 0 };

    double h=(alt/1000)*REARTH/((alt/1000)+REARTH);     //  geometric to geopotential altitude

    i=0; j=NTAB-1;  // starting values for binary search
    do
    {
        k=(i+j)/2;
        if (h < htab[k]) j=k; else i=k;
    }  while (j > i+1);

    double tgrad=gtab[i];                      // temp. gradient of local layer
    double tbase=ttab[i];                      // base temp. of local layer
    double deltah=h-htab[i];                   // height above local base
    double tlocal=tbase+tgrad*deltah;          // local temperature
    theta=tlocal/ttab[0];                                  // temperature ratio

    if (0.0 == tgrad)                                         // pressure ratio
        delta=ptab[i]*exp(-GMR*deltah/tbase);
    else
        delta=ptab[i]*pow(tbase/tlocal, GMR/tgrad);

    sigma=delta/theta;                                        //  density ratio
}

US76::~US76(){
}

double US76::p() {
    return PZERO*delta;
}

double US76::T() {
    return TZERO*theta;
}

double US76::rho(){
    return RHOZERO*sigma;
}