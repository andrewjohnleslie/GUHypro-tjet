/*!
 *  \brief      US76 atmosphere model
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

//This is an excerpt (modified) from PDAS, the original code can be found
//online at http://www.pdas.com/atmosdownload.html

#include <cmath>
#include <iostream>

class US76 {
private:
    const double TZERO          = 288.15;      ///< sea level temperature, kelvins
    const double PZERO          = 101325.0;        ///< sea-level pressure, N/sq.m
    const double RHOZERO        = 1.225;           ///< sea level density, kg/cu.m

    double delta;
    double theta;
    double sigma;

public:
    ///Constructor
    /**
     * @param alt altitude
     */
    US76(double alt);
    virtual ~US76();

    double p(); ///<Returns pressure
    double T(); ///<Returns temperature
    double rho(); ///<Returns density
};
