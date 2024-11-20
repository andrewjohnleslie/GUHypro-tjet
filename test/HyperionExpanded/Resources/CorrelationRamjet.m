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
close all
clear all

foot =  0.3048; %conversion feet to meters
lbm = 0.45359; %conversion libre to kg
g0 = 9.81;
inch = 0.0254;
rankine = 1/1.8;

%% Create an atmosphere model

% load geometry file 
vehicle = spaceplane('fullVehicle.stl');

% define some geometry parameters
vehicle.referenceArea = 300;  % [m^2]
vehicle.referenceLength = 20; % [m]
vehicle.wingspan = 16;        % [m]

% define the upper limit altitude of the model
boundaryAltitude = 300e3; % [m]

% create atmosphere model interpolant from US76 model
atmos = atmosphere(boundaryAltitude,'US76',vehicle);

%% Comparison
Alt = 56000*foot;
p = atmos.Query(Alt,'pressure');
T = atmos.Query(Alt,'temperature');
M = 3.5;
        
[Thrust, mfr, Data] = Hyperion(1.05*p,1.02*T,M,1);
Isp = Thrust/(sum(mfr)*g0);
gas = HyPro.mixture;
nodes = HyPro.Node.NodesInit(gas,Data);
nodes = nodes([2,5,10,11]);
            
SCCREAM.A = [27	8.24	22.5	44.3]*foot^2;
SCCREAM.M = [3.02 0.63 0.90 2.08];
SCCREAM.U = [3278 1083 2488 4870]*foot;
SCCREAM.p0 = [100.4 78.2 45.3 45.3]*lbm*g0/inch^2;
SCCREAM.T0 = [1345 1345 3446 3446]*rankine;
SCCREAM.Isp = 2370 + (Alt/foot - 5e4)*(2430 - 2370)/(7e4 - 5e4);

fprintf(1,'\n%5s','A');
fprintf(1,'%10g',SCCREAM.A);
fprintf(1,'\n%5s',' ');
fprintf(1,'%10g',nodes.A);

fprintf(1,'\n%5s','M');
fprintf(1,'%10g',SCCREAM.M);
fprintf(1,'\n%5s',' ');
fprintf(1,'%10g',nodes.M);

fprintf(1,'\n%5s','U');
fprintf(1,'%10g',SCCREAM.U);
fprintf(1,'\n%5s',' ');
fprintf(1,'%10g',nodes.U);

fprintf(1,'\n%5s','p0');
fprintf(1,'%10g',SCCREAM.p0);
fprintf(1,'\n%5s',' ');
fprintf(1,'%10g',nodes.p0);

fprintf(1,'\n%5s','T0');
fprintf(1,'%10g',SCCREAM.T0);
fprintf(1,'\n%5s',' ');
fprintf(1,'%10g',nodes.T0);

fprintf(1,'\n%5s','Isp');
fprintf(1,'%10g',SCCREAM.Isp);
fprintf(1,'\n%5s',' ');
fprintf(1,'%10g',Isp);
fprintf(1,'\n');
