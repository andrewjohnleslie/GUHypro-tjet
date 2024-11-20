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
AltFeet = [0,20,40,60]*1e3;
Label = cell(size(AltFeet));
for i=1:length(AltFeet)
    Label{i} = sprintf('Altitude=%.0f ft',AltFeet(i));
end
Altitude = AltFeet*foot;
M = [0.01,0.5,1,1.5,2,2.5,2.9];
Thrust = nan(length(M),length(Altitude));
Isp = Thrust;
mfr = Thrust;
Ct = Thrust;
Data = cell(length(M),length(Altitude));
gas = HyPro.mixture;

%multiSystem.selected = 3;
for j=1:length(Altitude)
    for i=1:length(M)
        p = atmos.Query(Altitude(j),'pressure');
        T = atmos.Query(Altitude(j),'temperature');
        
        fprintf(1,'H=%g M=%g\n',Altitude(j),M(i));
        tic
        try
            [Thrust(i,j), mfrX, Data{i,j}] = Hyperion(p,T,M(i),0);
            mfr(i,j) = sum(mfrX);
            N = HyPro.Node.NodesInit(gas,Data{i,j});
            Ct(i,j) = Thrust(i,j)/(0.5*N(1).rho*N(1).U^2*N(2).A);
        catch err
            if (strcmp(err.identifier,'HYPRO:error'))
                warning('HYPRO:error',err.message);

                Thrust(i,j) = NaN;
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
SCCREAM.Isp=[...
    409.3 409.3 445.0 561.9 828.3 1150.0 1335;
    380.0 389.8 425.5 510.0 633.4 828.3 977.7;
    367.1 367.1 389.8 445.0 519.7 600.9 698.4;
    354.1 354.1 367.1 383.3 415.8 464.5 513.2];

SCCREAM.Thrust = [...
    93.57 93.57 103.0 133.3 213.6 346.9 475.4;
    86.01 86.01 95.46 118.1 150.3 213.6 279.8;
    79.40 80.34 86.01 98.30 115.3 140.8 172.0;
    77.50 77.50 79.40 83.18 92.63 103.0 119.1];

SCCREAM.Ct = SCCREAM.Thrust.*Ct'./Thrust';

figure('Name','Thrust', 'NumberTitle','off', 'WindowStyle','docked')

A = imread('Hyperion_Olds1997','png');
X0 = [168,573]; X1 = [1051,44];
A = A(X0(2):-1:X1(2),X0(1):X1(1),:);
image([0,3],[0,500],A)
h(1) = gca;
set(gca,'YDir','normal','NextPlot','add',...
    'LineStyleOrder',{'-s','-^','-o','-*'})
SCCREAM.h(:,1) = plot(h(1),M,SCCREAM.Thrust,'k','DisplayName',Label);
legend(SCCREAM.h(:,1));

figure('Name','Spec. Thrust', 'NumberTitle','off', 'WindowStyle','docked')
set(gca,'NextPlot','add',...
    'LineStyleOrder',{'-s','-^','-o','-*'})
h(2) = gca;
SCCREAM.h(:,2) = plot(gca,M,SCCREAM.Ct,'k','DisplayName',Label);
legend(SCCREAM.h(:,2));

figure('Name','Specific Impulse', 'NumberTitle','off', 'WindowStyle','docked')

A = imread('HyperionIsp_Olds1997','png');
image(A)
X0 = [122,466]; X1 = [814,35];
A = A(X0(2):-1:X1(2),X0(1):X1(1),:);
image([0,3],[0,1400],A)
h(3) = gca;
set(gca,'YDir','normal','NextPlot','add',...
    'LineStyleOrder',{'-s','-^','-o','-*'})
SCCREAM.h(:,3) = plot(h(3),M,SCCREAM.Isp,'k','DisplayName',Label);
legend(SCCREAM.h(:,3));

%% plot data
plot(h(1),M,Thrust,'r','DisplayName',Label)
plot(h(2),M,Ct,'r','DisplayName',Label)
plot(h(3),M,Isp,'r','DisplayName',Label)