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
%GPstatLoop

clear all
close all

Ngen = 40;
for i=0:Ngen
    name = ['Generation-',num2str(i),'.dat'];
    data(i+1) = GPstat(name);
    
    best(i+1) = min(data(i+1).ft);
    average(i+1) = mean(data(i+1).ft(~isinf(data(i+1).ft)));
    errors(i+1) = 100*sum(data(i+1).status==2)/length(data(i+1).status);
    negative(i+1) = 100*sum(data(i+1).status==1)/length(data(i+1).status);
end

plot(0:Ngen,best,0:Ngen,average)

figure
plot(0:Ngen,errors,0:Ngen,negative)