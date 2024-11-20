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
%GP stat
% clear all
% close all

function data = GPstat(name,pieChart)

if nargin < 2
    pieChart = false;
end

fid = fopen(name);

C = textscan(fid,'%d%*[^:]%*1c%[^\n]%*1c%f%[^\n]%*1c');

fclose(fid);

filter = '-'; % 'Inlet';
id = strfind(C{2},filter);
for i=1:length(id)
    I(i) = ~isempty(id{i});
end
C{1} = C{1}(I);
C{2} = C{2}(I);
C{3} = C{3}(I);
C{4} = C{4}(I);

data.ID = C{1};
data.struct = C{2};
data.ft = C{3};
data.message = C{4};
data.status = zeros(size(C{1})); % 0 OK, 1 negative thrust, 2 error

TF(:,1) = strcmp('',C{4});
label{1} = 'Working';

msg = 'Error: flow choked in module:';
TF(:,2) = strncmp(msg,C{4},length(msg));
label{2} = 'Missing Feedback';
data.status(TF(:,2)) = 2;

msg = 'The thrust or the mfr is negative.';
TF(:,3) = strcmp(msg,C{4});
label{3} = 'Negative thrust';
data.status(TF(:,3)) = 1;

msg = 'Error: reduce method undefined.';
TF(:,4) = strncmp(msg,C{4},length(msg));
label{4} = 'Wrong feedback';
data.status(TF(:,4)) = 2;

msg = 'Error: X component cannot be negative.';
TF(:,5) = strncmp(msg,C{4},length(msg));
label{5} = 'Wrong X';
data.status(TF(:,5)) = 2;

msg = 'Error: function evaluations';
TF(:,6) = strncmp(msg,C{4},length(msg));
label{6} = 'Numerical';
data.status(TF(:,6)) = 2;

msg = 'Subsonic non adapted exhaust';
TF(:,7) = strncmp(msg,C{4},length(msg));
label{7} = 'Subsonic exhaust';
data.status(TF(:,7)) = 2;

msg = 'Error: flow always choked in module:';
TF(:,8) = strncmp(msg,C{4},length(msg));
label{8} = 'Always choked';
data.status(TF(:,8)) = 2;

msg = 'Error: A = ';
TF(:,9) = strncmp(msg,C{4},length(msg));
label{9} = 'A out of range';
data.status(TF(:,9)) = 2;

TF(:,10) = ~any(TF,2);
label{10} = 'Others';
data.status(TF(:,10)) = 2;

if pieChart
    q = sum(TF,1);
    
    tot = sum(q);
    if tot~=length(C{4})
        warning('Warning: total count is %d over %d elements',tot,length(C{4}));
    end
    
    pie3(q,label)
end
