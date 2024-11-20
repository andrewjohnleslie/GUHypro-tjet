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
Altitude = [0, 6e4*foot];
M = [0, 2.5];

%% Initialize figures
LineStyleOrder = {'-o','-*','-s','-v'};
figure('Name','Pressure', 'NumberTitle','off', 'WindowStyle','docked')
h(1) = axes('NextPlot','add','YScale','log');
ColorOrder = get(h(1),'ColorOrder');

figure('Name','Mach', 'NumberTitle','off', 'WindowStyle','docked')
h(2) = axes('NextPlot','add');

figure('Name','Area', 'NumberTitle','off', 'WindowStyle','docked')
h(3) = axes('NextPlot','add');

k = 1;
for i=1:length(Altitude)
    for j=1:length(M)
        Alt = round(Altitude(i)/100)*100;
        name = sprintf('h=%.0fm M=%.1f',Alt,M(j));
        p = atmos.Query(Altitude(i),'pressure');
        T = atmos.Query(Altitude(i),'temperature');
        try
            [Thrust, mfr, Data] = Hyperion(p,T,M(j),0);
            gas = HyPro.mixture;
            nodes = HyPro.Node.NodesInit(gas,Data);
            nodes = [nodes(1:6),nodes(8:end)];
        catch err
            if (strcmp(err.identifier,'HYPRO:error'))
                warning('h=%f M=%f \n%s',Altitude(i),M(j),err.message)
                Thrust = [];
                mfr = [];
                gas = mixture;
                nodes = HyPro.Node(gas);
            else
                rethrow(err);
            end
        end
        
        %% plot data
        plot(h(1),[nodes.p],LineStyleOrder{k},'DisplayName',name,...
            'Color',ColorOrder(k,:),'MarkerSize',15)
        
        plot(h(2),[nodes.M],LineStyleOrder{k},'DisplayName',name,...
            'Color',ColorOrder(k,:),'MarkerSize',15);
        
        % figure('Name','MFR', 'NumberTitle','off', 'WindowStyle','docked')
        % h(3) = gca;
        % for i=1:length(nodeList)
        %     mfr(i) = nodeList(i).mfr;
        % end
        % plot(mfr)
        
        % figure('Name','Composition', 'NumberTitle','off', 'WindowStyle','docked')
        % h(4) = gca;
        % plot(reshape([nodeList.X],4,length(nodeList))')
        
        plot(h(3),[nodes.A],LineStyleOrder{k},'DisplayName',name,...
            'Color',ColorOrder(k,:),'MarkerSize',15)
        plot(h(3),[nodes.Amax],'p','DisplayName',name,...
            'Color',ColorOrder(k,:),'MarkerSize',15)
        
        % figure('Name','Velocity', 'NumberTitle','off', 'WindowStyle','docked')
        % h(6) = axes('NextPlot','add');
        % for i=1:length(nodeList)
        %     U(i) = nodeList(i).U;
        % end
        % plot(U)
        
        % figure('Name','Total Pressure', 'NumberTitle','off', 'WindowStyle','docked')
        % h(7) = axes('NextPlot','add');
        % for i=1:length(nodeList)
        %     p0(i) = nodeList(i).p0;
        % end
        % plot(p0)
        
        k = k + 1;
    end
end

for i=1:length(h)
    DH = 0.12;
    set(h(i),'XTickLabel',[],'XTick',1:length(nodes),'XGrid','on',...
        'OuterPosition',[0,DH,1,1-DH],'Box','on');
    
    Yl = get(h(i),'YLim');
    t = text(1:length(nodes),Yl(1)*ones(1,length(nodes)),nodes.Name,...
        'HorizontalAlignment','right','VerticalAlignment','top', ...
        'Rotation',45,'Parent',h(i),'FontSize',20);
end

%% Comparison
% n = length(inlet.nodes);
% n = 2;
% Node = [2,n+1,n+5,n+6];
% Mach = [3.02, 0.63, 0.9, 2.08];
% U = [3278, 1083, 2488, 4870]*foot;
% p0 = [100.4, 78.2, 45.3, 45.3]*6894;
% plot(h(2),Node,Mach,'or')
% plot(h(6),Node,U,'or')
% plot(h(7),Node,p0,'or')