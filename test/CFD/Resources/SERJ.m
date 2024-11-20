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

foot = 0.3048; %0.34; %conversion feet to meters
lbm = 0.45359; %conversion libre to kg
g0 = 9.81;
inch = 0.0254;
psi = 6894.75728; %conversion psi to Pa

%% Thermal model
load('/home/trb12187/Matlab/SpaceplaneModel/Fluidynamic/ThermoModels/PredefinedModels/janafModels.mat');
janafModels = janafModels([1,3,2,4]);

gas = mixture('Air+Fuel',janafModels);

%% Inputs
freeStream = node();
freeStream.p = 101325;
freeStream.T = 300;
freeStream.gas = gas;
freeStream.X('O2','N2') = [0.209,1-0.209];   %Air
freeStream.M = 0.75;
freeStream.A = NaN;

freeStreamMod = node();
freeStreamMod.gas = freeStream.gas;
exit = node();

%% Forebody
forebody = wedge('forebody',freeStream,freeStreamMod);
forebody.delta = 8;
forebody.calculate;

%% System
multiSystem = multiMode('SERJ',freeStreamMod,exit);

system = multiSystem.add('Supercharged Ejector');
%Nodes 1   2     3      4            5      6      7    8     9    10     11     12    13
A =   [70, 70,   52,   43.49-5.24, 43.39, 43.49, 43.49, 0, 5.24, 76.24, 76.24, 76.24, 107]*foot^2;
%Nodes      1            2        3        4           5         6     7
name = {'Pre Intake','Intake','Throat','Pinch Point','Fan','Primary','Mixer End',...
    'Rocket Throat','Rocket Outlet','End Diffuser','Injection','End Chamber','Nozzle'};
    %        8             9              10           11             12       13
system.initNodes(length(A));
for i=1:length(A)
    system.nodes(i).A = A(i);
    system.nodes(i).gas = gas;
    system.nodes(i).name = name{i};
end

%% Inlet
inlet = system.add(propModules.adaptedInlet,'inlet',[1,2,3,4]);
inlet.chokingFeedback = inlet;

%% Fan
system.add(propModules.isentropicDuct,'post pinch point',[4,5]);

compress = system.add(propModules.compressor,'Fan',[5,6]);
compress.R = 1.3;

%% Primary - ejector
primaryNozzle = system.add(propModules.isentropicDuct,'primary',[8,9]);
primaryNozzle.choked = true;
mfr = 392*lbm;
p0 = 1500*psi;
T0 = 3675; %adiabatic flame calculation
primaryNozzle.N1.M = 1;
primaryNozzle.N1.T = 1000; % Only to asses gamma
primaryNozzle.N1.X('H2O') = 1;
[pstr_p0,Tstr_T0] = isentropic(primaryNozzle.N1.M,primaryNozzle.N1.gamma);
primaryNozzle.N1.p = p0*pstr_p0;
primaryNozzle.N1.T = T0*Tstr_T0;
primaryNozzle.N1.A = mfr/(primaryNozzle.N1.rho*primaryNozzle.N1.a);

primary = system.add(propModules.mixer,'primary',[6,7,9]);
primary.teta = 0;
primary.verbose = 1;
primary.chokingFeedback = inlet; %preInlet;
primary.eff = 1;

%% Secondary - Post combustor
system.add(propModules.isentropicDuct,'diffuser',[7,10]);

secondary = system.add(propModules.injectionPhi,'secondary',[10,11]);
secondary.Phi = 1;
secondary.Xf.gas = secondary.N1.gas;
secondary.Xf('H2') = 1;
secondary.chokingFeedback = secondary;

combustor = system.add(propModules.combChamber,'combustor',[11,12],1,1);
combustor.Cf = 0;
combustor.chokingFeedback = secondary;

%% Nozzle
nozzle = system.add(propModules.adaptedNozzle,'nozzle',[12,13]);
nozzle.choked = true;
nozzle.Aemax = nozzle.N2.A;

% in case of inlet spilling (especially in supersonic) better to adapt the
% nozzle to the inlet conditions.
nozzle.freeStream = inlet.N1;

%% Pure Ramjet mode
ramjet = multiSystem.add('Ramjet mode',1);

inlet1 = ramjet.exchange(1,propModules.adaptedInlet,'Inlet',[1,2,3,4]);
% inlet1.nodes(3).A = 27*foot^2;
% inlet1.nodes(3).gas = gas;
inlet1.chokingFeedback = inlet1;

primary1 = ramjet.exchange(3,propModules.neutralLink,'Primary',[5,6]);

%% Execution
multiSystem.selected = 1;
tic
multiSystem.calculate;
toc

%% plot data
Thrust = system.thrust/(1e3*lbm*g0)

% nodeList = [freeStream, multiSystem.modes{multiSystem.selected}.modules{1}.nodes, system.nodes([3:4,6:end])];
nodeList = [freeStream, system.nodes([1:6,8:end])];

figure('Name','Pressure', 'NumberTitle','off', 'WindowStyle','docked')
plot([nodeList.p])
h(1) = gca;

figure('Name','Mach', 'NumberTitle','off', 'WindowStyle','docked')
h(2) = axes('NextPlot','add');
plot([nodeList.M]);

figure('Name','MFR', 'NumberTitle','off', 'WindowStyle','docked')
h(3) = gca;
for i=1:length(nodeList)
    mfr(i) = nodeList(i).mfr;
end
plot(mfr)

figure('Name','Composition', 'NumberTitle','off', 'WindowStyle','docked')
h(4) = gca;
plot(reshape([nodeList.X],4,length(nodeList))')

figure('Name','Area', 'NumberTitle','off', 'WindowStyle','docked')
h(5) = gca;
plot([nodeList.A])

figure('Name','Velocity', 'NumberTitle','off', 'WindowStyle','docked')
h(6) = axes('NextPlot','add');
for i=1:length(nodeList)
    U(i) = nodeList(i).U;
end
plot(U)

figure('Name','Total Pressure', 'NumberTitle','off', 'WindowStyle','docked')
h(7) = axes('NextPlot','add');
for i=1:length(nodeList)
    p0(i) = nodeList(i).p0;
end
plot(p0)

for i=1:length(h)
    set(h(i),'XTickLabel',{nodeList.name},'XTick',1:length(nodeList));
end

%% Comparison
% n = length(inlet.nodes);
n = 2;
Node = [2,n+1,n+5,n+6];
Mach = [3.02, 0.63, 0.9, 2.08];
U = [3278, 1083, 2488, 4870]*foot;
p0 = [100.4, 78.2, 45.3, 45.3]*6894;
plot(h(2),Node,Mach,'or')
plot(h(6),Node,U,'or')
plot(h(7),Node,p0,'or')