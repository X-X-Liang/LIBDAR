function zprime = Kirkwood_Bethe_method(t,z,y,varargin)
% This function tracks the outgoing characteristics in the liquid around
% a spherical bubble using the Kirkwood-Bethe hypothesis see 
% Section 3.3. Acoustic and shock wave emission in the Ref: 
% 
% Liang et al. J. Fluid. Mech. 940, A5 (2022). DOI: 10.1017/jfm.2022.202
%
%  Usage:
%  zprime = Kirkwood_Bethe_method(t,z,y)
%
%  Input arguments:
%  t    - time axis for calculating the propagation of characteristics
%  z(1) - u, local liquid velocity
%  z(2) - r, radial distance
%  y    - y = r(h+u^2/2), quantity that propagates outward
%
%  Optional parameter
%  p0     - ambient pressure % Update version 1.1
%
%  Output arguments:
%  zprime(1) - du/dt, acceleration of local liquid
%  zprime(2) - dr/dt = u, local liquid velocity

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
%% Parse optional input arguments
p = inputParser;
defaultp0     = 1e5;    % default value for ambient pressure 1e5 Pa
addParameter(p,'p0',defaultp0);
parse(p,varargin{:});   % Parse param.-value pairs
param = p.Results;      % transfer res. to structure
clear p                 % delete the object p that is not needed
p0     = param.p0;      % ambient pressure as an input parameter

% mass density and speed of sound at ambient pressure p0
[rho, c0] = TaitEOS(p0); 
%% parameters for the calculations:
% p0    = 1e5;          % ambient pressure
% rho   = 998;          % density of water 998kg/m3
% c0    = 1483;         % speed of sound in water 1483m/s

B     = 3.14e8;       % 314MPa, constant in Tait EOS
n     = 7;            % constant in Tait EOS
p_inf = p0;           % pressure in the liquid, far away from the bubble; here assumed to be equal to p0

p = (p_inf + B)*((y/z(2) - z(1)^2/2)*(n-1)*rho/(n*(p_inf + B))+1)^(n/(n-1))-B; % Eq.(3.33) in Liang2022
c = c0*((p + B)/(p_inf+B))^((n-1)/2/n);                                        % Eq.(3.32) in Liang2022

zprime(1) = 1/(c-z(1))*((c+z(1))*y/z(2)^2-2*c^2*z(1)/z(2));                    % Eq.(3.31) in Liang2022
zprime(2) = z(1)+c;                                                                

zprime = zprime(:);