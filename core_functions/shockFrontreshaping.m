function [x_output,y_output,xout,index_R,index_x0] = shockFrontreshaping(x,y)
% This function reshapes the multi-valued u(r) and p(r) curves at shock front. See Ref:
%
% Liang et al. J. Fluid. Mech. 940, A5 (2022). DOI: 10.1017/jfm.2022.202
%
% Usage:
% [x_output,y_output,xout,index_R,index_x0] = shockFrontreshaping(x,y)
%
% input arguments
% x - radius r
% y - velocity u, or pressure p
%
% output arguments
% x_output - radius r after reshaping
% y_output - u or p after reshaping
%
%   LICENSE DISCLAIMER:
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program. If not, see <http://www.gnu.org/licenses/>.
%
%   Copyright (C) 2024  Xiao-Xuan(Joe) Liang
%   Insitute of Biomedical Optics, University of Luebeck, 
%   Peter-Monnik-Weg 4, 23562 Luebeck, Germany
%   Correspondence to x.liang@uni-luebeck.de

%% Get rid of NaNs
x(isnan(x)) = 0;            % we must get rid of the NaN's
y(isnan(y)) = 0;            % we must get rid of the NaN's
j = find(x);                % and only find the valid values
x = x(j);
y = y(j);


%% Find Right point 
index_R = find(diff(x)<0,1,'first');  % simplest way to find the R point


%% Find Left point - (x0,y0) 
%  the turning point in the f2 and f3 curves, this must be right, otherwise it's the whole programme would be problematic

[~,index_x0] = min(x(index_R:end)); % simplest way to find the L point in L2&L3 curves
index_x0     = index_x0+index_R-1;  % update index_x0 to the whole curve
x0           = x(index_x0);
y0           = y(index_x0);

%% Delineate line f1
index_left_f1 = find(x(1:index_R)<x0,1,'last');
if isempty(index_left_f1) % make sure there is an index for the left_f1 point 11.03.2025
    index_left_f1 = 1;
end
x_line_f1     = x(index_left_f1:index_R);
y_line_f1     = y(index_left_f1:index_R);
f1            = griddedInterpolant(x_line_f1,y_line_f1,'pchip');

%% Delineate line f2
xx  = x(index_R:end);
yy  = y(index_R:end);
t0  = yy > y0;
t00 = xx > x0;
t1  = t0 & t00; % line f2
x_line_f2=[x0;flipud(xx(t1))];
y_line_f2=[y0;flipud(yy(t1))];


%  Bad points in line f2: Zigzag points in line f2 make the line non-monotonical, thus shall be dropped
while (find(diff(x_line_f2)<0))  % judge if it is monotonical by using the x-array;
    badpoints_line_f2            = find(diff(x_line_f2)<0);
    x_line_f2(badpoints_line_f2) = [];
    y_line_f2(badpoints_line_f2) = [];
end

f2 = griddedInterpolant(x_line_f2, y_line_f2, 'pchip');

%% Delineate line f3
if x(end)<x(index_R) % for expansion phase
    x_line_f3 = [x0;(x0+x(index_R))/2;x(index_R)];
    y_line_f3 = [y0;y0;y0];
else
    t2 = ~t0 & t00;
    x_line_f3=[x0;xx(t2)];
    y_line_f3=[y0;yy(t2)];
    %  Bad points in line f3: Zigzag points in line f2 make the line non-monotonical, thus shall be dropped
    while (find(diff(x_line_f3)<0))
        badpoints_line_f3=find(diff(x_line_f3)<0);
        x_line_f3(badpoints_line_f3) = [];
        y_line_f3(badpoints_line_f3) = [];
    end
end   
f3 = griddedInterpolant(x_line_f3, y_line_f3, 'pchip');

%% Integration f1, f2 and f3, find xout for equal area
xout = fzero(@(X)integral2(@(x,y)ones(size(x)),X,x(index_R),@(x)f2(x),@(x)f1(x))...
               - integral2(@(x,y)ones(size(x)),x0,X,@(x)f3(x),@(x)f2(x)),[x0,x(index_R)]);

%% Output reshaped curve           
index_x1_out = find(x(1:index_R)<xout,1,'last');
y1_out = f1(xout);
index_x2_out = find(flipud(x)<xout,1,'first')-1;  % turn the x array up side down
index_x2_out = length(x)-index_x2_out+1;
y2_out = f3(xout);

if x(end)<x(index_R)
    x_output=[x(1:index_x1_out); xout;   xout;   xout];
    y_output=[y(1:index_x1_out); y1_out; y2_out; y(end)];
else
    x_output=[x(1:index_x1_out); xout;   xout;   x(index_x2_out:end)];
    y_output=[y(1:index_x1_out); y1_out; y2_out; y(index_x2_out:end)];
end

end