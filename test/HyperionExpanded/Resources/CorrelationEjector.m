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
Alt = 0;
p = atmos.Query(Alt,'pressure');
T = atmos.Query(Alt,'temperature');
M = 0.5;
        
[Thrust, mfr, Data] = Hyperion(p,T,M,0);
gas = HyPro.mixture;
nodes = HyPro.Node.NodesInit(gas,Data);
nodes = nodes([3,5,11,12]);
            
SCCREAM.A = [27	8.24	22.5	17.3]*foot^2;
SCCREAM.M = [0.5	0.57 NaN		1.51];
SCCREAM.U = [558.1	636.5 NaN		5290.1]*foot;
SCCREAM.p0 = [17.7	17.7	50.8	50.8]*lbm*g0/inch^2;
SCCREAM.T0 = [544.6	544.6	5544.7	5544.7]*rankine;

fprintf(1,'\nA ');
fprintf(1,'%10g',SCCREAM.A);
fprintf(1,'\n  ');
fprintf(1,'%10g',nodes.A);

fprintf(1,'\nM ');
fprintf(1,'%10g',SCCREAM.M);
fprintf(1,'\n  ');
fprintf(1,'%10g',nodes.M);

fprintf(1,'\nU ');
fprintf(1,'%10g',SCCREAM.U);
fprintf(1,'\n  ');
fprintf(1,'%10g',nodes.U);

fprintf(1,'\np0 ');
fprintf(1,'%10g',SCCREAM.p0);
fprintf(1,'\n   ');
fprintf(1,'%10g',nodes.p0);

fprintf(1,'\nT0 ');
fprintf(1,'%10g',SCCREAM.T0);
fprintf(1,'\n   ');
fprintf(1,'%10g',nodes.T0);
fprintf(1,'\n');

x = 1:length(nodes.M);
figure
plot(x,[SCCREAM.M;nodes.M])
figure
plot(x,[SCCREAM.U;nodes.U])
figure
plot(x,[SCCREAM.p0;nodes.p0])
figure
plot(x,[SCCREAM.T0;nodes.T0])
