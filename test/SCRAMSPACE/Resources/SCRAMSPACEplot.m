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

gas = HyPro.mixture;
out = SCRAMSPACE(0.4);
N = HyPro.Node.NodesInit(gas,out);

out = SCRAMSPACE(0);
Nfroz = HyPro.Node.NodesInit(gas,out);

%% read data
CFD.react = dlmread('LorrainCFD.csv',',',[1,0,282,1]);
CFD.froz = dlmread('LorrainCFD.csv',',',[284,0,482,1]);

%% plot data
figure('Name','wall Pressure', 'NumberTitle','off', 'WindowStyle','docked')

xCFD = 0.17937;

plot(CFD.react(:,1),1000*CFD.react(:,2),'r','DisplayName','CFD reacting')
hold on
plot(CFD.froz(:,1),1000*CFD.froz(:,2),'b','DisplayName','CFD frozen')

x = [-0.179,0,0,0.15,0.38,0.755-0.179]+xCFD;
plot(x,N.p,'r','DisplayName','HyPro reacting')
plot(x,Nfroz.p,'b','DisplayName','HyPro frozen')

fid = fopen('/home/trb12187/OpenFOAM/trb12187-2.1.1/run/Propulsion/Lorrain_Roconstr/NonComb-1/0.0003858049/oneDim.dat');
fscanf(fid,'%[^\n]%1c',2);
A = fscanf(fid,'%f',[6,inf]);
fclose(fid);

plot(A(1,:),A(2,:),'--b','DisplayName','CFD average')

figure('Name','Speed', 'NumberTitle','off', 'WindowStyle','docked')
plot(x,N.U)

hold on

plot(A(1,:)-xCFD,A(5,:))

figure('Name','Temperature', 'NumberTitle','off', 'WindowStyle','docked')
plot(x,N.T)

hold on

plot(A(1,:)-xCFD,A(3,:))