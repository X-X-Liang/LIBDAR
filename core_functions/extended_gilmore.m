function xprime = extended_gilmore(t,x,R0,Rnbd,Rnc1,Rnc2,tRmax,tRmax2,tp,varargin)
% This function used the extended Gilmore model to calculate the temporal
% development of the bubble radius and pressure inside the bubble, 
% which considers the rapid increase of bubble wall velocity during 
% laser-induced energy deposition in a compressible liquid.
% see Section 3.1 and 3.2. in the Ref: 

% Liang et al. J. Fluid. Mech. 940, A5 (2022). DOI: 10.1017/jfm.2022.202

% and Section II A. in the Ref:
% Liang&Vogel,https://doi.org/10.48550/arXiv.2501.13749
%
% Usage:
% xprime = extended_gilmore(t,x,R0,Rnbd,Rnc1,Rnc2,tRmax,tRmax2,tp)
%
% Input arguments:
% t      - time axis for calculating bubble dynamics
% x(1)   - R bubble radius 
% x(2)   - U velocity of bubble wall 
% R0     - initial bubble radius
% Rnbd   - equilibrium bubble radius after optical breakdown
% Rnc1   - equilibrium bubble radius at the first bubble collapse
% Rnc2   - equilibrium bubble radius at the second bubble collapse
% tRmax  - time instant at Rmax
% tRmax2 - time instatn at Rmax2
% tp     - pulse duration
%
% Optional parameter
% RNP    - radius of extra hard core, e.g. nanoparticles, default =0 
% tRmax3 - time instant at Rmax3, default = infinite
% Twall  - temperature at the bubble wall, default = 20oC
% 'JSC'  - choice of the jump-start condition (JSC). Options:
%                'no' - no jump-start condition
%                'first-order'- first-order approximation for inertial
%                confined condition
%                'second-order'- second-order approximation for inertial
%                confined condition
%                'general' - second-order approximation for general
%                condition (default)
% p0     - ambient pressure % Update version 1.1
%
% Output arguments:
% xprime(1) - time derivatives of R, bubble wall velocity
% xprime(2) - time derivatives of U, bubble wall acceleration
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

%% Parse optional input arguments
p = inputParser;
defaultJSC = 'general'; % default value for jump-start-condition (JSC)
defaultRNP = 0;         % default value for the radius of hard core, e.g. nanoparticles
defaulttRmax3 = 1000;   % default value for tRmax, set a large value as "infinite", in seconds
defaultTwall  = 293;    % default value for the temperature at bubble wall,in Kelvin
defaultRnc3   = Rnc2;   % default value for Rncoll3
defaultisIsotherm  = 1; % default value for the kappa value for late oscillations t > tRmax3, 1 for isothermal, 0 for adiabatic
defaultp0     = 1e5;    % default value for ambient pressure 1e5 Pa

validationFcn = @(str) any(strcmp(str,{'no','first-order','second-order','general'}));
addParameter(p,'JSC',defaultJSC,validationFcn);
addParameter(p,'RNP',defaultRNP);
addParameter(p,'tRmax3',defaulttRmax3);
addParameter(p,'Twall', defaultTwall);
addParameter(p,'Rnc3', defaultRnc3);
addParameter(p,'isIsotherm', defaultisIsotherm);
addParameter(p,'p0',defaultp0);

parse(p,varargin{:});   % Parse param.-value pairs
param = p.Results;      % transfer res. to structure
clear p                 % delete the object p that is not needed

Twall  = param.Twall;
tRmax3 = param.tRmax3;
RNP    = param.RNP;
Rnc3   = param.Rnc3;
p0     = param.p0;      % ambient pressure as an input parameter

% mass density and speed of sound at ambient pressure p0
[rho, c0] = TaitEOS(p0);  
%% constants
% p0      = 1e5;          % ambient pressure 1e5 Pa at STP
% c0      = 1483;         % speed of sound 1483m/s at STP
% rho     = 998;          % density of water 998kg/m3 at STP
kappa   = 4/3;          % adiabatic exponent 4/3 for H2O vapor 
B       = 3.14e8;       % parameter in Tait EOS
A       = B+p0;         % B+p0 parameters for the Tait EOS
n       = 7;            % default 7
c1      = 5190;         % constants in the Hugoniot curve see Eq.(3.11) in Liang2022
c2      = 25306;        % constants in the Hugoniot curve see Eq.(3.11) in Liang2022
p_inf   = p0;           % pressure in the liquid, far away from the bubble; here assumed to be equal to p0
sigma   = 0.072583;     % surface tension of water in N/m 0.072583 at room temperature
mu      = 1.046e-3;     % viscosity in 1mPa s =1.046e-3 N/m2 s for water at room temperature

%% Temperature-dependent surface tension and viscosity for shock heated bubble wall, see Fig. 11 in Liang2022JFM
if t > tRmax3     % switch at t = t_Rmax3, kappa = 1 for linear oscillations
    sigma = sigma_T(Twall); % T-dependent surface tension
    mu    = mu_T(Twall);    % T-dependent viscosity
    if param.isIsotherm == 1
        kappa = 1;          % isothermal condition
    else
        % adiabatic
    end
end

%% Temperoal change of the equilibrium bubble radius Rn(t)
Rnt = (R0^3 + (Rnbd^3-R0^3)/2/tp*(t-tp/pi*sin(pi/tp*t)))^(1/3);  % Rnt during pulse duration
Rnt = Rnt.*(t<=2*tp) + Rnbd*(t>2*tp & t<=tRmax) +...
      Rnc1*(t>tRmax & t<=tRmax2) + Rnc2*(t>tRmax2 & t<=tRmax3)+Rnc3*(t>tRmax3);              % Rnt switch at tRmax1 and tRmax2

%% Consideration of van der Waals hard core Rvdw and other hard cores, e.g. gold nanoparticles RNP
if RNP > 0             % extra hard core like NPs are considered
    Rvdw = RNP*(t<tRmax)+(((1/9)^3*(Rnt.^3-RNP^3)+RNP^3).^(1/3)).*(t>=tRmax); % old version Rvdw = RNP,now RNP+shell of vdw gas
    assert(R0 > RNP);  % radius of the initial hot vapor size must be larger than the nanoparticle size.
else                   % no extra hard core
    Rvdw = 0*(t<tRmax) + Rnt/9*(t>=tRmax);
end

%% Quantities at bubble wall
pn = p0 + 2*sigma/Rnt;   
p  = pn*((Rnt^3-Rvdw^3)/(x(1)^3-Rvdw^3))^kappa -...
     2*sigma/x(1) -4*mu*x(2)/x(1);      % pressure at the bubble wall Eq. (3.3) in Liang2022JFM
H  = n/(n-1)*A^(1/n)/rho*((p+B)^((n-1)/n)-(p_inf+B)^((n-1)/n));  % liquid enthalpy Eq. (3.6) in Liang2022JFM
dHdR = 1/rho*(A/(p+B))^(1/n)*...        % Eq. (3.7) in Liang2022JFM
    (-3*pn*kappa*x(1)^2*(Rnt^3-Rvdw^3)^(kappa)/((x(1)^3-Rvdw^3)^(kappa+1))+2*sigma/x(1)^2+4*mu*x(2)/x(1)^2);   
C = (c0^2+(n-1)*H)^0.5;                 % speed of sound at the bubble wall Eq. (3.5) in Liang2022JFM

%% Jump start condition
if t<=2*tp % jump start condition applies to the time window of energy deposition 
    switch param.JSC
        case 'no' % no jump-start
            % disp('no jump-start condition is used')
            dUpdt = 0;
        case 'first-order'  % first-order approximation for inertial confined condition, Eq. 3.23 in Liang2022JFM
            % disp('first-order approximation for inertial confined condition is used')
            dUpdt = (2*p0*Rnt+3*sigma)*(Rnbd^3-R0^3)*(1-cos(pi*t/tp))/(3*R0^4*tp*rho*c0);
        case 'second-order' % 2nd order approximation for inertial confined condition,   Eq. 3.30 in Liang2022JFM
            % disp('second-order approximation for inertial confined condition is used')
            ps_bd = p0*Rnt^4/R0^4+2*sigma*Rnt^3/R0^4;                                  % Eq. 3.19 in Liang2022JFM
            dUpdt = (2*p0*Rnt+3*sigma)*(Rnbd^3-R0^3)*(1-cos(pi*t/tp))/(3*R0^4*tp)/sqrt(rho^2*c0^2+4*rho*c2*ps_bd/(log(10)*c1));
        case 'general'      % general case
            % disp('general condition is used')
            p_gas  = pn*((Rnt^3-Rvdw^3)/(x(1)^3-Rvdw^3))^kappa;
            Rntdot = (Rnbd^3-R0^3)*(1-cos(pi*t/tp))/6/tp/Rnt^2;
            Pdot  = 3*kappa*(p0+2*sigma/Rnt)*(Rnt^3-Rvdw^3)^kappa/(x(1)^3-Rvdw^3)^kappa*(Rnt^2*Rntdot/(Rnt^3-Rvdw^3)-x(1)^2*x(2)/(x(1)^3-Rvdw^3))-...
                2*sigma*Rntdot/Rnt^2*(Rnt^3-Rvdw^3)^kappa/(x(1)^3-Rvdw^3)^kappa; % +2*sigma/x(1)^2+4*mu*(x(2)/x(1)^2-xprime(2)/x(1)); surface tension and viscosity
            %   Pdot   = 4*p0*Rnt^3/x(1)^4*(Rntdot-x(2)*Rnt/x(1))+2*sigma*Rnt^2/x(1)^4*(3*Rntdot-4*Rnt*x(2)/x(1)); % for water RNP = 0 and kappa = 4/3
            if Pdot < 0     % shock wave pressure should accelerate not deaccelerate bubble movement
                Pdot = 0;
            end
            dUpdt = Pdot/sqrt(rho^2*c0^2+4*rho*c2*p_gas/(log(10)*c1));
    end
else
    dUpdt = 0;
end

%% Rate equations
xprime(1) = x(2);                                        % speed of bubble wall U=dR/dt  
xprime(2) = -3/2*x(2)^2/x(1)*(C-x(2)/3)/(C-x(2))+...
            H/x(1)*(C+x(2))/(C-x(2))+x(2)/C*dHdR+dUpdt;  % Gilmore equation
xprime    = xprime(:);

end