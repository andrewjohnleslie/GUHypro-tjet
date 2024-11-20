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
STX(1).Name = 'STX-05';
STX(2).Name = 'STX-06';
STX(3).Name = 'STX-09';
STX(1).Number = 5;
STX(2).Number = 6;
STX(3).Number = 9;
STX(1).Color = 'b'; 
STX(2).Color = 'g';
STX(3).Color = 'm';
STX(1).Vt = dlmread('Traj_V_time.csv',',',[1,0,115,1]);
STX(2).Vt = dlmread('Traj_V_time.csv',',',[117,0,222,1]);
STX(3).Vt = dlmread('Traj_V_time.csv',',',[224,0,286,1]);
STX(1).ZS = 1e3*dlmread('Traj_Z_S.csv',',',[1,0,117,1]);
STX(2).ZS = 1e3*dlmread('Traj_Z_S.csv',',',[119,0,234,1]);
STX(3).ZS = 1e3*dlmread('Traj_Z_S.csv',',',[236,0,289,1]);
STX(1).ZS30 = 1e3*dlmread('Traj_Z_S.csv',',',[291,0,291,1]);
STX(2).ZS30 = 1e3*dlmread('Traj_Z_S.csv',',',[293,0,293,1]);
STX(3).ZS30 = 1e3*dlmread('Traj_Z_S.csv',',',[295,0,295,1]);
STX(1).ZM = dlmread('Traj_Z_M.csv',',',[1,0,65,1]);
STX(2).ZM = dlmread('Traj_Z_M.csv',',',[67,0,117,1]);
STX(3).ZM = dlmread('Traj_Z_M.csv',',',[119,0,153,1]);
for i=1:length(STX)
    STX(i).ZM(:,2) = 1e3*STX(i).ZM(:,2);
    [~,I] = sort(STX(i).ZM(:,2));
    STX(i).ZM = STX(i).ZM(I,:);
end

STX(1).Data.pi = dlmread('InternalPress.csv',',',[61,0,113,1]);
STX(2).Data.pi = dlmread('InternalPress.csv',',',[115,0,170,1]);
STX(3).Data.pi = dlmread('InternalPress.csv',',',[172,0,209,1]);
for i=1:length(STX)
    STX(i).Data.pi(:,2) = 1e5*STX(i).Data.pi(:,2);
end

%% Reinterpolate trajectory
TrajMach = true;

for i=1:length(STX)
    STX(i).ZS(1,:) = [0,0];
    
    STX(i).ZS30(2) = interp1(STX(i).ZS(:,1),STX(i).ZS(:,2),STX(i).ZS30(1));
    if all(STX(i).ZS(:,1)~=STX(i).ZS30(1))
        I = STX(i).ZS(:,1)>STX(i).ZS30(1);
        STX(i).ZS = [STX(i).ZS(~I,:); STX(i).ZS30; STX(i).ZS(I,:)];
    end
    
    L = cumsum(sqrt(diff(STX(i).ZS(:,1)).^2 + diff(STX(i).ZS(:,2)).^2));
    STX(i).ZS = [STX(i).ZS, [0; L]];
    STX(i).ZS30(3) = STX(i).ZS(STX(i).ZS(:,1)==STX(i).ZS30(1),3);
    L = cumtrapz(STX(i).Vt(:,1),STX(i).Vt(:,2));
    L30 = interp1(STX(i).Vt(:,1),L,30);
    L = L - L30 + STX(i).ZS30(3);
    STX(i).Vt = [STX(i).Vt, L];
    
    Z = interp1(STX(i).ZS(:,3),STX(i).ZS(:,2),STX(i).Vt(:,3));
    STX(i).Vt = [STX(i).Vt, Z];
    
    t = 20:10:90;
    Z = interp1(STX(i).Vt(:,1),STX(i).Vt(:,4),t);
    fprintf(1,'%10.0f ',t); fprintf(1,'\n');
    fprintf(1,'%10.0f ',Z); fprintf(1,'\n');
    fprintf(1,'\n');
    
    STX(i).t = STX(i).Vt(:,1);
    STX(i).V = STX(i).Vt(:,2);
    STX(i).Alt = STX(i).Vt(:,4);
    if TrajMach
        STX(i).Mach = interp1(STX(i).ZM(:,2),STX(i).ZM(:,1),STX(i).Alt);
        I = ~isnan(STX(i).Mach);
        STX(i).t = STX(i).t(I);
        STX(i).V = STX(i).V(I);
        STX(i).Alt = STX(i).Alt(I);
        STX(i).Mach = STX(i).Mach(I);
    end
end

%% Atmosphere correction
TempCorr = false;

gas = HyPro.mixture('Test/STATALTEXRamjet',{'JETA','O2','N2','H2O','CO2'});
freeStream = HyPro.Node(gas);
freeStream.X([0,0.209,1-0.209,0,0]);

for i=1:length(STX)
    for j=1:length(STX(i).Alt)
        STX(i).p(j) = atmos.Query(STX(i).Alt(j),'pressure');
        STX(i).T(j) = atmos.Query(STX(i).Alt(j),'temperature');
        freeStream.p(STX(i).p(j));
        freeStream.T(STX(i).T(j));
        
        if TempCorr
            if ~TrajMach
                error('Error: TempCorr can be selected only if TrajMach is selected');
            end
            a = STX(i).V(j)/STX(i).Mach(j);
            STX(i).T(j) = fzero(@(T)invertA(T,freeStream,a),freeStream.T);
        end

        if TrajMach && ~TempCorr
            STX(i).V(j) = STX(i).Mach(j)*freeStream.a;
        else
            freeStream.U(STX(i).V(j));
            STX(i).Mach(j) = freeStream.M;
        end
    end
    
    if TempCorr
        STX(i).p(1) = atmos.Query(STX(i).Alt(1),'pressure');
        for j=1:length(STX(i).Alt)-1
            freeStream.p(STX(i).p(j));
            freeStream.T(STX(i).T(j));
            STX(i).p(j+1) = STX(i).p(j) - freeStream.rho*g0*(STX(i).Alt(j+1) - STX(i).Alt(j));
        end
    end
end

%% Calculate Hyperion with HyPro
pCorr = false; %Activate free stream pressure correction by means of internal pressure
Aref = pi*(0.263/2)^2;
Iint = 5;

% STX(i).p(j) = STX(i).p(j)*0.8;

for i=1:length(STX)
    for j=1:length(STX(i).Alt)
        p = STX(i).p(j); 
        T = STX(i).T(j); 
        
        try
            [STX(i).Thrust(j), mfrX, Data] = STATALTEXRamjet(p,T,STX(i).Mach(j),STX(i).Number,0);
            STX(i).mfr(j) = sum(mfrX);
            N = HyPro.Node.NodesInit(gas,Data);
            STX(i).zita(j) = N(2).mfr/(N(1).rho*N(1).U*Aref);
            STX(i).rho(j) = N(1).rho;
            STX(i).pi(j) = N(Iint).p;
            
            if pCorr 
                pref = interp1(STX(i).Data.pi(:,1),STX(i).Data.pi(:,2),STX(i).t(j));
                for k=1:5
                    G = 1;
                    p = p*G*pref/STX(i).pi(j);
                    
                    [STX(i).Thrust(j), mfrX, Data] = STATALTEXRamjet(p,T,STX(i).Mach(j),STX(i).Number,0);
                    STX(i).mfr(j) = sum(mfrX);
                    N = HyPro.Node.NodesInit(gas,Data);
                    STX(i).zita(j) = N(2).mfr/(N(1).rho*N(1).U*Aref);
                    STX(i).rho(j) = N(1).rho;
                    STX(i).pi(j) = N(Iint).p;
                    
                    if (pref/STX(i).pi(j) <= 1.01) && (pref/STX(i).pi(j) >= 0.99)
                        break
                    end
                    if k==5
                        warning('Warning: Maximum number of iteration reached')
                    end
                end
                STX(i).p(j) = p;
            end
            
            STX(i).Thrust(j) = 0.96*N(end).A*(N(end).rho*N(end).U^2 + N(end).p) - N(1).A*(N(1).rho*N(1).U^2 + N(1).p) ...
                - N(1).p*(N(end).A - N(1).A);
%             STX(i).Thrust(j) = STX(i).Thrust(j) - N(end).A*(N(end).rho*N(end).U^2 + N(end).p) + 0.96*N(end).A*(N(end).rho*N(end).U^2 + N(end).p);
            STX(i).Ct(j) = STX(i).Thrust(j)/(0.5*N(1).rho*STX(i).V(j)^2*Aref);
        catch err
            if (strcmp(err.identifier,'HYPRO:error'))
                warning('HYPRO:error',err.message);
                
                STX(i).Thrust(j) = NaN;
                STX(i).mfr(j) = NaN;
                STX(i).zita(j) = NaN;
                STX(i).rho(j) = NaN;
                STX(i).Ct(j) = NaN;
                STX(i).pi(j) = NaN;
            else
                rethrow(err);
            end
        end
        
        STX(i).Isp(j) = STX(i).Thrust(j)/(STX(i).mfr(j)*g0);
        %     fprintf(1,'Thrust=%f Isp=%f\n',STX(i).Thrust(j),STX(i).Isp(j));
        %     fprintf(1,'\n');
    end
end

%% Read data
STX(1).Data.Isp = dlmread('Isp.csv',',',[1,0,9,1]);
STX(2).Data.Isp = dlmread('Isp.csv',',',[11,0,18,1]);
STX(3).Data.Isp = dlmread('Isp.csv',',',[20,0,35,1]);
IspUp = dlmread('Isp.csv',',',[55,0,132,1]);
IspDown = dlmread('Isp.csv',',',[134,0,208,1]);

STX(1).Data.Thrust = dlmread('Thrust.csv',',',[1,0,21,1]);
STX(2).Data.Thrust = dlmread('Thrust.csv',',',[23,0,42,1]);
STX(3).Data.Thrust = dlmread('Thrust.csv',',',[44,0,66,1]);
for i=1:length(STX)
    STX(i).Data.Thrust(:,2) = 10*STX(i).Data.Thrust(:,2);
end

STX(1).Data.zita = dlmread('InletProperties.csv',',',[1,0,46,1]);
STX(2).Data.zita = STX(1).Data.zita;
STX(3).Data.zita = STX(1).Data.zita;

%% Plot
figure('Name','Trajectory', 'NumberTitle','off', 'WindowStyle','docked')
subplot(2,2,1,'NextPlot','add','box','on')
for i=1:length(STX)
    plot(STX(i).Mach,STX(i).Alt,'Color',STX(i).Color,'DisplayName',['HyPro ',STX(i).Name]);
    plot(STX(i).ZM(:,1),STX(i).ZM(:,2),'+','Color',STX(i).Color,'DisplayName',STX(i).Name)
end
subplot(2,2,2,'NextPlot','add','box','on')
Alt = linspace(0,6e4,100);
for i=1:length(Alt)
    p(i) = atmos.Query(Alt(i),'pressure');
    T(i) = atmos.Query(Alt(i),'temperature');
end
plot(T,Alt,'k--','DisplayName','STD Atmosphere')
for i=1:length(STX)
    plot(STX(i).T,STX(i).Alt,'Color',STX(i).Color,'DisplayName',['HyPro ',STX(i).Name]);
end
subplot(2,2,3,'NextPlot','add','box','on','XScale','log')
plot(p,Alt,'k--','DisplayName','STD Atmosphere')
for i=1:length(STX)
    plot(STX(i).p,STX(i).Alt,'Color',STX(i).Color,'DisplayName',['HyPro ',STX(i).Name]);
end

figure('Name','Isp', 'NumberTitle','off', 'WindowStyle','docked')
axes('NextPlot','add','box','on');
for i=1:length(STX)
    h(i) = plot(STX(i).Mach,STX(i).Isp,'Color',STX(i).Color,'DisplayName',['HyPro ',STX(i).Name]);
    plot(STX(i).Data.Isp(:,1),STX(i).Data.Isp(:,2),'+','Color',STX(i).Color,'DisplayName',STX(i).Name)
end
plot(IspDown(:,1),IspDown(:,2),'--','DisplayName','Variation')
plot(IspUp(:,1),IspUp(:,2),'--','DisplayName','Variation')

figure('Name','Isp Error', 'NumberTitle','off', 'WindowStyle','docked')
axes('NextPlot','add','box','on');
for i=1:length(STX)
    I = find(diff(STX(i).Mach)>=0,1,'last');
    Isp = interp1(STX(i).Mach(1:I),STX(i).Isp(1:I),STX(i).Data.Isp(:,1));
    plot(STX(i).Data.Isp(:,1),100*abs(STX(i).Data.Isp(:,2)-Isp)./STX(i).Data.Isp(:,2),'+','Color',STX(i).Color,'DisplayName',STX(i).Name)
end

figure('Name','Thrust', 'NumberTitle','off', 'WindowStyle','docked')
axes('NextPlot','add','box','on');
for i=1:length(STX)
    plot(STX(i).t,STX(i).Thrust,'Color',STX(i).Color,'DisplayName',['HyPro ',STX(i).Name]);
%     plot(STX(i).t,STX(i).Thrust1,'--','Color',STX(i).Color,'DisplayName',['HyPro ',STX(i).Name]);
    plot(STX(i).Data.Thrust(:,1),STX(i).Data.Thrust(:,2),'+','Color',STX(i).Color,'DisplayName',STX(i).Name)
end

figure('Name','Thrust Error', 'NumberTitle','off', 'WindowStyle','docked')
h1 = axes('NextPlot','add','box','on');
figure('Name','Thrust Error', 'NumberTitle','off', 'WindowStyle','docked')
h2 = axes('NextPlot','add','box','on');
for i=1:length(STX)
    T = interp1(STX(i).t,STX(i).Thrust,STX(i).Data.Thrust(:,1));
    Mach = interp1(STX(i).t,STX(i).Mach,STX(i).Data.Thrust(:,1));
    Alt = interp1(STX(i).t,STX(i).Alt,STX(i).Data.Thrust(:,1));
    plot(h1,Mach,100*abs(STX(i).Data.Thrust(:,2)-T)./STX(i).Data.Thrust(:,2),'+','Color',STX(i).Color,'DisplayName',STX(i).Name)
    plot(h2,Alt,100*abs(STX(i).Data.Thrust(:,2)-T)./STX(i).Data.Thrust(:,2),'+','Color',STX(i).Color,'DisplayName',STX(i).Name)
end

figure('Name','Thrust', 'NumberTitle','off', 'WindowStyle','docked')
axes('NextPlot','add','box','on');
for i=1:length(STX)
    plot(STX(i).Mach,STX(i).Ct,'Color',STX(i).Color,'DisplayName',['HyPro ',STX(i).Name]);
    dyn = STX(i).Thrust./STX(i).Ct;
    dyn = interp1(STX(i).t,dyn,STX(i).Data.Thrust(:,1));
    Mach = interp1(STX(i).t,STX(i).Mach,STX(i).Data.Thrust(:,1));
%     rho = interp1(STX(i).t,STX(i).rho,STX(i).Data.Thrust(:,1));
%     V = interp1(STX(i).t,STX(i).V,STX(i).Data.Thrust(:,1));
%     Ct = STX(i).Data.Thrust(:,2)./(0.5*Aref*rho.*V.^2);
    Ct = STX(i).Data.Thrust(:,2)./dyn;
    plot(Mach,Ct,'+','Color',STX(i).Color,'DisplayName',STX(i).Name)
end

figure('Name','Intake', 'NumberTitle','off', 'WindowStyle','docked')
axes('NextPlot','add','box','on');
for i=1:length(STX)
    plot(STX(i).Mach,STX(i).zita,'Color',STX(i).Color,'DisplayName',['HyPro ',STX(i).Name]);
    plot(STX(i).Data.zita(:,1),STX(i).Data.zita(:,2),'+','Color',STX(i).Color,'DisplayName',STX(i).Name)
end

figure('Name','Pressure', 'NumberTitle','off', 'WindowStyle','docked')
axes('NextPlot','add','box','on');
for i=1:length(STX)
    plot(STX(i).t,STX(i).pi,'Color',STX(i).Color,'DisplayName',['HyPro ',STX(i).Name]);
    plot(STX(i).Data.pi(:,1),STX(i).Data.pi(:,2),'+','Color',STX(i).Color,'DisplayName',STX(i).Name)
end

%% save trajectory data
fid = fopen('Traj_processed.csv','w');
for i=1:length(STX)
    fprintf(fid,'Name = %s\n',STX(i).Name);
    fprintf(fid,'Altitude,Pressure,Temperature,Mach\n');
    fprintf(fid,'%f,%f,%f,%f\n',[STX(i).Alt,STX(i).p',STX(i).T',STX(i).Mach]');
end
fclose(fid);