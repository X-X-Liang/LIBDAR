function [r_reshaped,u_reshaped,p_reshaped,rout,pout] = shockFronttransfer(r,u,p)
% This function is used to transfer the shock front position at reshape
% u(r) curves to the multi-valued p(r) curves. See Ref:
%
% Liang et al. J. Fluid. Mech. 940, A5 (2022). DOI: 10.1017/jfm.2022.202
%
% Usage:
% [r_reshaped,u_reshaped,p_reshaped,rout,pout] = shockFronttransfer(r,u,p)
%
% input arguments
% r - radius 
% u - velocity
% p - pressure
%
% output argument
% r_reshaped,u_reshaped,p_reshaped are reshaped curves
%(rout,pout) is the point at sw front that is used to determine the shock wave decay slope
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
r(isnan(r)) = 0;            % we must get rid of the NaN's
u(isnan(u)) = 0;            % we must get rid of the NaN's
p(isnan(p)) = 0;
j = find(r);                % and only find the valid values
r = r(j);
u = u(j);
p = p(j);

%% reshape u(r) curve
[r_reshaped,u_reshaped,xout,index_R,index_x0] = shockFrontreshaping(r,u);

%% reshaped p(r) curve according to the reshaped u(r) curve at xout 
% find wave front at f1 line
index_x1_out = find(r(1:index_R)<xout,1,'last');
y1_out = interp1(r(index_x1_out:index_x1_out+1),p(index_x1_out:index_x1_out+1),xout,'linear');
% find wave front at f3 line
if r(end)<r(index_R)   % for expansion phase to polish the shape of p(r)
    p_L = p(index_x0);
    y2_out = p_L;
    
else
    index_x2_out = find(flipud(r)<xout,1,'first')-1;  % turn the x array up side down
    index_x2_out = length(r)-index_x2_out+1;
    y2_out = interp1(r(index_x2_out-1:index_x2_out),p(index_x2_out-1:index_x2_out),xout,'linear');
end

if r(end)<r(index_R)
    p_reshaped=[p(1:index_x1_out); y1_out; y2_out; p(end)];
else
    p_reshaped=[p(1:index_x1_out); y1_out; y2_out; p(index_x2_out:end)];
end
rout = xout;
pout = y1_out;


end

