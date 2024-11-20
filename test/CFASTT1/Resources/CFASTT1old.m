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

%% Thermal model
load('/home/trb12187/Matlab/SpaceplaneModel/Fluidynamic/ThermoModels/PredefinedModels/janafModels.mat');
janafModels = janafModels([1,3,2,4]);

gas = mixture('Air+Fuel',janafModels);

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

%% Inputs
Altitude = 0; % 6e4*foot;
freeStream = node();
freeStream.p = atmos.Query(Altitude,'pressure'); %33.52;
freeStream.T = atmos.Query(Altitude,'temperature'); %255.7;
freeStream.gas = gas;
freeStream.X('O2','N2') = [0.209,1-0.209];   %Air
freeStream.M = 3.5;
freeStream.A = NaN;

freeStreamMod = node();
freeStreamMod.gas = freeStream.gas;
exit = node();

%% System
multiSystem = multiModeMachP0('Hyperion ERJ',freeStream,exit);

system = multiSystem.add('Ejector Mode');
%Nodes 1   2     3      4      5      6     7           8     9    10    11    12
A =   [27, 27, 0.25*27, 8.24, 11.25, 11.25, NaN, 11.25-8.24, 22.5, 22.5, 22.5, 95]*foot^2;
name = {'Pre Intake','Intake','Throat','Pinch Point','Primary','Mixer End',...
    'Dummy','Rocket Outlet','End Diffuser','Injection','End Chamber','Nozzle'};
system.initNodes(length(A));
for i=1:length(A)
    system.nodes(i).A = A(i);
    system.nodes(i).gas = gas;
    system.nodes(i).name = name{i};
end
system.nodes(2).Amax = A(2);
system.nodes(3).Amax = A(3);
system.nodes(3).Amin = 0;
system.nodes(12).Amax = A(12);

%% Inlet
inlet = system.add(propModules.adaptedThroatInlet,'inlet',[1,2,3,4]);
inlet.chokingFeedback = inlet;

%% Primary - ejector
system.add(propModules.isentropicDuct,'post pinch point',[4,5]);

rocket = system.add(propModules.rocket,'Rocket',[7,8]);
rocket.N2.X('H2O') = 1;
rocket.N2.M = 6;
rocket.mfr = 216*lbm;
rocket.Isp = 462;

primary = system.add(propModules.mixer,'primary',[5,6,8]);
primary.verbose = 1;
primary.chokingFeedback = inlet; %preInlet;
primary.eff = 0.67;

%% Secondary - Post combustor
system.add(propModules.isentropicDuct,'diffuser',[6,9]);

secondary = system.add(propModules.injectionPhi,'secondary',[9,10]);
secondary.PhiMax = 1;
secondary.Xf.gas = secondary.N1.gas;
secondary.Xf('H2') = 1;
secondary.chokingFeedback = secondary;

combustor = system.add(propModules.combChamber,'combustor',[10,11],1,1);
combustor.Cf = 0;
combustor.chokingFeedback = secondary;

%% Nozzle
nozzle = system.add(propModules.effAdaptedNozzle,'nozzle',[11,12]);
nozzle.choked = true;
nozzle.p0Ratio(1);

% in case of inlet spilling (especially in supersonic) better to adapt the
% nozzle to the inlet conditions.
nozzle.freeStream = inlet.N1;

%% Pure Ramjet mode
ramjet = multiSystem.add('Ramjet mode',1);

ramjet.nodes(3) = ramjet.nodes(3).clone; %detach node #3 from previous mode
inlet1 = ramjet.exchange(1,propModules.adaptedInlet,'Inlet',[1,2,3,4]);
inlet1.Nt.Amax = inlet1.Nin.Amax;
inlet1.Nt.A = inlet1.Nt.Amax;
inlet1.chokingFeedback = inlet1;

primary1 = ramjet.exchange(4,propModules.neutralLink,'primary',[5,6]);
ramjet.modules = [ramjet.modules(1:2),ramjet.modules(4:end)];

%% Pure Rocket mode
rocketMode = multiSystem.add('Rocket mode',1);

primary2 = rocketMode.exchange(4,propModules.isentropicDuct,'primary',[8,6]);

rocketMode.exchange(6,propModules.neutralLink,'secondary',[9,11]);

nozzle2 = rocketMode.exchange(8,propModules.effAdaptedNozzle,'nozzle',[11,12]);
nozzle2.choked = false;

rocketMode.modules = [rocketMode.modules(3:6),rocketMode.modules(8)];
rocketMode.nodes(1) = rocketMode.nodes(1).clone; %detach node #1 from previous mode
rocketMode.nodes(1).A = 0; 
rocketMode.N1 = rocketMode.nodes(1);
nozzle2.freeStream = rocketMode.N1;

rocketMode.nodes = rocketMode.nodes([1,6:9,11:end]);

%% Off mode
off = multiSystem.add('closed engine');
off.initNodes(2);
off.add(propModules.neutralLink,'dummy',[1,2]);

%% set ranges and save
multiSystem.rangesM = [  [0,3];     [2.5,8];   [0,inf]; [0,inf] ];
multiSystem.rangesP0 = [ [1e4,inf]; [5e5,inf]; [0,inf]; [0,inf] ];

% engine = multiSystem;
% save propCFASTT1 engine

%% Execution
tic
multiSystem.calculate;
toc

%% plot data
Thrust = multiSystem.thrust
Isp = multiSystem.specificImpulse
[PropMfr,PropMfrX] = multiSystem.propellantMfr

switch multiSystem.autoSelect.moduleName
    case 'Ejector Mode'
        nodeList = [freeStream, system.nodes([1:6,9:end])];
    case 'Ramjet mode'
        nodeList = [freeStream, multiSystem.autoSelect.nodes([1:6,9:end])];
    case 'Rocket mode'
        nodeList = [freeStream, multiSystem.autoSelect.nodes([1:2,4:end])];
    otherwise
        nodeList = [];
end

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
h(5) = axes('NextPlot','add');
plot([nodeList.A])
plot([nodeList.Amax],'vr')

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