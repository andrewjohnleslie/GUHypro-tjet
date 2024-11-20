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
clear all 
close all

foot =  0.3048; %conversion feet to meters
lbm = 0.45359; %conversion libre to kg
g0 = 9.81;
inch = 0.0254;

%% Create an atmosphere model

% load geometry file 
fuselage = spaceplane('fuselageWing.stl');

% define some geometry parameters
fuselage.referenceArea = 300;  % [m^2]
fuselage.referenceLength = 20; % [m]
fuselage.wingspan = 16;        % [m]

% define the upper limit altitude of the model
boundaryAltitude = 200e3; % [m]

% create atmosphere model interpolant from US76 model
atmos = atmosphere(boundaryAltitude,'US76',fuselage);

%% Read Data
A = csvread('Hyperion_traj.csv',1,0);
M = A(:,1);
Altitude = A(:,2)*foot;

SCCREAM.Isp = dlmread('SCCREAM_verif_Isp.csv',',',[1,0,30,1]);
Astronautics.Isp = dlmread('SCCREAM_verif_Isp.csv',',',[32,0,58,1]);
MarquardtERJ.Isp = dlmread('SCCREAM_verif_Isp.csv',',',[60,0,63,1]);
MarquardtESJ.Isp = dlmread('SCCREAM_verif_Isp.csv',',',[64,0,68,1]);

SCCREAM.Ct = dlmread('SCCREAM_verif_Ct.csv',',',[1,0,54,1]);
SCCREAM.M = M;
SCCREAM.CtM = interp1(SCCREAM.Ct(:,1),SCCREAM.Ct(:,2),M);
Astronautics.Ct = dlmread('SCCREAM_verif_Ct.csv',',',[107,0,136,1]);
Astronautics.M = M;
Astronautics.CtM = interp1(Astronautics.Ct(:,1),Astronautics.Ct(:,2),M);
MarquardtERJ.Ct = dlmread('SCCREAM_verif_Ct.csv',',',[56,0,77,1]);
MarquardtESJ.Ct = dlmread('SCCREAM_verif_Ct.csv',',',[78,0,105,1]);

%% Calculate Hyperion with HyPro
gas = HyPro.mixture;
Aref = 27.0*foot^2;
for j=1:length(Altitude)
    p = atmos.Query(Altitude(j),'pressure');
    T = atmos.Query(Altitude(j),'temperature');
    
    if M(j)>2.8
        mode = 1;
    else
        mode = 0;
    end
    
    fprintf(1,'H=%g M=%g\n',Altitude(j),M(j));
    tic
    try
        [Thrust(j), mfrX, Data{j}] = Hyperion(p,T,M(j),mode);
        mfr(j) = sum(mfrX);
        N = HyPro.Node.NodesInit(gas,Data{j});
        q(j) = 0.5*N(1).rho*N(1).U^2;
        Ct(j) = Thrust(j)/(q(j)*Aref);
    catch err
        if (strcmp(err.identifier,'HYPRO:error'))
            warning('HYPRO:error',err.message);
            
            Thrust(j) = NaN;
            mfr(j) = NaN;
        else
            rethrow(err);
        end
    end
    toc
    
    Isp(j) = Thrust(j)/(mfr(j)*g0);
    Thrust(j) = Thrust(j)/(1e3*lbm*g0);
    fprintf(1,'Thrust=%f Isp=%f\n',Thrust(j),Isp(j));
    fprintf(1,'\n');
end

%% Plot
figure('Name','Isp', 'NumberTitle','off', 'WindowStyle','docked')
h(1) = plot(SCCREAM.Isp(:,1),SCCREAM.Isp(:,2),'k','DisplayName','SCCREAM');
hold on
h(2) = plot(Astronautics.Isp(:,1),Astronautics.Isp(:,2),'DisplayName','Astronautics');
h(3) = plot(Astronautics.Isp(1:2:end,1),Astronautics.Isp(1:2:end,2),'sq','DisplayName','Astronautics');
h(4) = plot(NaN,NaN,'-sq','DisplayName','Astronautics');
h(5) = plot(MarquardtERJ.Isp(:,1),MarquardtERJ.Isp(:,2),'-o','DisplayName','Marq. ERJ');
h(6) = plot(MarquardtESJ.Isp(:,1),MarquardtESJ.Isp(:,2),'-v','DisplayName','Marq. ESJ');
h(7) = plot(M,Isp,'--k','DisplayName','HyPro');
set(h,'LineWidth',1,'MarkerSize',8)
set(h(2:6),'Color',[0.5,0.5,0.5])
legend([h(1),h(4:end)])

figure('Name','Thrust', 'NumberTitle','off', 'WindowStyle','docked')
h(1) = plot(SCCREAM.Ct(:,1),SCCREAM.Ct(:,2),'k','DisplayName','SCCREAM');
hold on
h(2) = plot(Astronautics.Ct(:,1),Astronautics.Ct(:,2),'DisplayName','Astronautics');
h(3) = plot(Astronautics.Ct(1:2:end,1),Astronautics.Ct(1:2:end,2),'sq','DisplayName','Astronautics');
h(4) = plot(NaN,NaN,'-sq','DisplayName','Astronautics');
h(5) = plot(MarquardtERJ.Ct(:,1),MarquardtERJ.Ct(:,2),'DisplayName','Marq. ERJ');
h(6) = plot(MarquardtERJ.Ct(1:2:end,1),MarquardtERJ.Ct(1:2:end,2),'o','DisplayName','Marq. ERJ');
h(7) = plot(NaN,NaN,'-o','DisplayName','Marq. ERJ');
h(8) = plot(MarquardtESJ.Ct(:,1),MarquardtESJ.Ct(:,2),'DisplayName','Marq. ESJ');
h(9) = plot(MarquardtESJ.Ct(1:2:end,1),MarquardtESJ.Ct(1:2:end,2),'v','DisplayName','Marq. ESJ');
h(10) = plot(NaN,NaN,'-v','DisplayName','Marq. ESJ');
h(11) = plot(M(M>2.8),Ct(M>2.8),'--k','DisplayName','HyPro');
set(h,'LineWidth',1,'MarkerSize',8)
set(h(2:10),'Color',[0.5,0.5,0.5])
legend([h(1),h(4),h(7),h(10:end)])

figure('Name','Trajectory', 'NumberTitle','off', 'WindowStyle','docked')
plot(M,q*foot^2/(lbm*g0))

figure('Name','Trajectory', 'NumberTitle','off', 'WindowStyle','docked')
plot(M,Altitude)