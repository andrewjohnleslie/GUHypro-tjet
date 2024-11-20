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
    
%% Load the fuselage of Skylon

% load geometry file 
fuselage = spaceplane('fuselageWing.stl');

% define some geometry parameters
fuselage.referenceArea = 300;  % [m^2]
fuselage.referenceLength = 20; % [m]
fuselage.wingspan = 16;        % [m]

%% Create an atmosphere model
% define the upper limit altitude of the model
boundaryAltitude = 200e3; % [m]

% create atmosphere model interpolant from US76 model
atmos = atmosphere(boundaryAltitude,'US76',fuselage);

%% Execution
Aref = 27.0*foot^2;
gas = HyPro.mixture;
AltFeet = [30,50,70,110]*1e3;
Label = cell(size(AltFeet));
for i=1:length(AltFeet)
    Label{i} = sprintf('Altitude=%.0f ft',AltFeet(i));
end
Altitude = AltFeet*foot;
M = 2:0.5:6;
Thrust = nan(length(M),length(Altitude));
Isp = Thrust;
Ct = Thrust;
mfr = Thrust;
Data = cell(length(M),length(Altitude));

%multiSystem.selected = 3;
for j=1:length(Altitude)
    for i=1:length(M)
        p = atmos.Query(Altitude(j),'pressure');
        T = atmos.Query(Altitude(j),'temperature');
        
        fprintf(1,'H=%g M=%g\n',Altitude(j),M(i));
        tic
        try
            [Thrust(i,j), mfrX, Data{i,j}] = Hyperion(p,T,M(i),1);
            nodes = HyPro.Node.NodesInit(gas,Data{i,j});
            Ct(i,j) = Thrust(i,j)/(0.5*nodes(1).rho*nodes(1).U^2*Aref);
            mfr(i,j) = sum(mfrX);
        catch err
            if (strcmp(err.identifier,'HYPRO:error'))
                warning('HYPRO:error',err.message);

                Thrust(i,j) = NaN;
                Ct(i,j) = NaN;
                mfr(i,j) = NaN;
            else
                rethrow(err);
            end
        end
        toc
        fprintf(1,'\n');
        
        Isp(i,j) = Thrust(i,j)/(mfr(i,j)*g0);
        Thrust(i,j) = Thrust(i,j)/(1e3*lbm*g0);
    end
end

%% Image read
SCCREAM.Isp = [3070 3270 2500 2200 2500 2600 2550 2500 2450;
               3107 3510 2700 2370 2620 2780 2740 2680 2630;
               3170 3730 2900 2430 2680 2840 2800 2750 2680;
               3020 3520 3150 2600 2800 2900 2850 2790 2730];
SCCREAM.Ct = [0.42 0.64 0.79 1.27 1.6 1.55 1.5 1.45 1.4;
              0.42 0.68 0.88 1.35 1.7 1.7 1.66 1.63 1.58;
              0.42 0.62 0.94 1.43 1.74 1.71 1.64 1.59 1.54;
              0.42 0.77 1 1.57 1.8 1.75 1.7 1.65 1.6];
SCCREAM.Thrust = SCCREAM.Ct.*Thrust'./Ct';


figure('Name','Thrust', 'NumberTitle','off', 'WindowStyle','docked')

% A = imread('Hyperion_Olds1997','png');
% X0 = [168,573]; X1 = [1051,44];
% A = A(X0(2):-1:X1(2),X0(1):X1(1),:);
% image([0,3],[0,500],A)
h(1) = gca;
set(gca,'YDir','normal','NextPlot','add',...
    'LineStyleOrder',{'-s','-^','-o','-*'})
SCCREAM.h(:,1) = plot(h(1),M,SCCREAM.Ct,'k','DisplayName',Label);
legend(SCCREAM.h(:,1));

figure('Name','Specific Impulse', 'NumberTitle','off', 'WindowStyle','docked')

% A = imread('HyperionIsp_Olds1997','png');
% image(A)
% X0 = [122,466]; X1 = [814,35];
% A = A(X0(2):-1:X1(2),X0(1):X1(1),:);
% image([0,3],[0,1400],A)
h(2) = gca;
set(gca,'YDir','normal','NextPlot','add',...
    'LineStyleOrder',{'-s','-^','-o','-*'})
SCCREAM.h(:,2) = plot(h(2),M,SCCREAM.Isp,'k','DisplayName',Label);
legend(SCCREAM.h(:,2));

figure('Name','Thrust', 'NumberTitle','off', 'WindowStyle','docked')
h(3) = gca;
set(gca,'YDir','normal','NextPlot','add',...
    'LineStyleOrder',{'-s','-^','-o','-*'})
SCCREAM.h(:,3) = plot(h(3),M,SCCREAM.Thrust,'k','DisplayName',Label);
legend(SCCREAM.h(:,3));

%% plot data
plot(h(1),M,Ct,'r','DisplayName',Label)
plot(h(2),M,Isp,'r','DisplayName',Label)
plot(h(3),M,Thrust,'r','DisplayName',Label)