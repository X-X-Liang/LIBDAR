function [tau_th,Rmax_th,Tosc_th, Results_Matrix] = inertial_conf_tau_loop(R0,Rnbd_R0_ratio,tau_Array,Inertial_threshold,varargin)
% This function calculate the change of bubble parameters as a function of laser pulse duration for a given R0. See Ref:
% Liang&Vogel2024 Preprint...
% Usage:
% [tau_th,Rmax_th,Tosc_th, Results_Matrix] =
% inertial_conf_tau_loop(R0,Rnbd_R0_ratio,tau_Array,Inertial_threshold)
%
% Explicit input arguments [in SI units]:
% R0 - initial bubble size
% Rnbd_R0_ratio - Rnbd/R0
% tau_Array - array of pulse duration
% Inertial_threshold - threshold for inertial confinement 
%
% Optional parameters
% RNP - radius of extra hard core,e.g. nanoparticles
% 
% Output arguments:
% tau_th  - pulse duration at inertial conf. threshold
% Rmax_th - Rmax at inertial conf. threshold
% Tosc_th - bubble oscillation time at inertial conf. threshold
% Results_Matrix - matrix to store results
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
%
close all;
%% Parse optional input arguments
p = inputParser;
defaultRNP    = 0;      % default value for the radius of hard core, e.g. nanoparticles
validScalarnonNegNum = @(x) isnumeric(x) && isscalar(x) && (x >= 0); % non-negative scalar number
addParameter(p,'RNP',defaultRNP,validScalarnonNegNum);
parse(p,varargin{:});   % Parse param.-value pairs
param = p.Results;      % transfer res. to structure
clear p                 % delete the object p that is not needed

RNP = param.RNP;   
%% Constants
p0    = 1e5;           % ambient pressure in Pa
sigma = 0.072583;      % surface tension 0.072583N/m
kappa = 4/3;           % adiabatic exponent ; polytropic exponent

% preallocate array
L_tau = length(tau_Array);
Results_Matrix = zeros(6,L_tau); % Row 1 Rmax, Row 2 Tosc, Row 3 R2tau/R0, Row 4 Pbd, Row 5 Ubd, Row 6 tau_L

t1   = 0;
t2   = 10*1.83*R0*2;   % estimate a sufficiently long calculation time
Rnbd = R0*Rnbd_R0_ratio;
Rnc  = Rnbd/4;         % For studying bubble expansion, the value of Rnc is not important. Here we use a value close to that in the Ref. Liang2022 JFM. In reality, Rnc depends on bubble size as well as plasma energy density.

for i=1:L_tau
    tau_L = tau_Array(i);

    % inital calculation assuming no condensation to find Rmax and tRmax                                                                   
    Rn  = Rnbd;
    Rn2 = Rnbd;
    Rn3 = Rnbd;
    tRmax = t2;
    tRmax2= t2;
   
    x0 = [R0 0];    % initial condition  
    opts = odeset('RelTol',1e-8,'AbsTol',1e-9,'NormControl','off'); 
    [t,x] = ode45(@(t,x)extended_gilmore(t,x,R0,Rn,Rn2,Rn3,tRmax,tRmax2,tau_L),[t1 t2],x0,opts);
    R = x(:,1);     % find Rmax 
    [~,mt] = max(R);
    t_Rmax = t(mt); % find tRmax

    % recalculating assuming Rnc=Rnbd/4 at tRmax
    Rn  = Rnbd;
    Rn2 = Rnc;
    Rn3 = Rnc;
    tRmax = t_Rmax;
    tRmax2= t2;
    [t,x] = ode45(@(t,x)extended_gilmore(t,x,R0,Rn,Rn2,Rn3,tRmax,tRmax2,tau_L),[t1 t2],x0,opts);
    R = x(:,1);
    U = x(:,2);

    % find R|t=2tau_L
    R_2tau = interp1(t,R,2*tau_L,"spline"); % query the value of R(t) at t=2tau_L.
    
    % find Rmax
    [Rmax,mt] = max(R);
   
    % find Tosc1
    [~,nt]    = min(R(mt:length(R)));% find the min
    index_Rmin= nt+mt-1; % index of min.

    [pks,~]   = findpeaks(R(1:index_Rmin),t(1:index_Rmin));
    if length(pks)==1 % to check if the min. is the first collapse
    else
        [~,locs2]      = findpeaks(R(1:index_Rmin));
        [~,index_Rmin] = min(R(locs2(1):locs2(2)));% find the minima inbetween the 1st peak and the 2nd peak
        index_Rmin     = index_Rmin+locs2(1)-1;
    end
    Tosc1     = t(index_Rmin); % 1st osc. time in second
    
    % find Ubd,max
    U_bd_max  = max(U(1:mt));  % max. expansion velocity
    
    % find Pbd,max
    Rnt = (R0^3 + (Rn^3-R0^3)/2/tau_L*(t-tau_L/pi*sin(pi/tau_L*t))).^(1/3); % the equilibrium pressure rises during the laser pulse from R0 to Rn and 
                                                                            % is switched at time tRmax to the value Rn2
    Rnt = Rnt.*(t<=2*tau_L) + Rn*(t>2*tau_L & t<tRmax) + Rn2*(t>tRmax & t<tRmax2)+Rn3*(t>tRmax2);
    
    if RNP > 0             % van der Waals hard core Rvdw considering RNP 2022.12.15
        Rvdw = RNP*(t<tRmax)+(((1/9)^3*(Rnt.^3-RNP^3)+RNP^3).^(1/3)).*(t>=tRmax); % before Rvdw = RNP
        assert(R0 > RNP);  % radius of the initial hot vapor size is larger than the nanoparticle size.
    else
        Rvdw=0*(t<tRmax)+Rnt/9.*(t>=tRmax);
    end
    
    pn = p0 + 2*sigma./Rnt;                              % gas pressure at the current value of the equilibrium radius, Rnt
    pg = pn.*((Rnt.^3-Rvdw.^3)./(R.^3-Rvdw.^3)).^kappa;  % gas pressure in the bubble at the bubble radius R = R(t) with van der Waals core
    P_bd_max = max(pg(1:mt));                            % max. expansion pressure in Pa.

% store results
    Results_Matrix(1,i) = Rmax; 
    Results_Matrix(2,i) = Tosc1;
    Results_Matrix(3,i) = R_2tau/R0;
    Results_Matrix(4,i) = P_bd_max;
    Results_Matrix(5,i) = U_bd_max;
    Results_Matrix(6,i) = tau_L;
end

%% find the threshold values
if max(Results_Matrix(3,:)) > Inertial_threshold 
    tau_th  = interp1(Results_Matrix(3,:), Results_Matrix(6,:), Inertial_threshold,'spline'); % Inertial_threshold R/R0 = 2^(1/3) 
    Tosc_th = interp1(Results_Matrix(6,:), Results_Matrix(2,:), tau_th,'spline'); 
    Rmax_th = interp1(Results_Matrix(6,:), Results_Matrix(1,:), tau_th,'spline'); 

    % Plots
    figure()
    subplot(5,1,1)
    semilogx(Results_Matrix(6,:),Results_Matrix(3,:),'g',tau_th,Inertial_threshold,'o','LineWidth',1.5);
    ylabel('R|_{t=2tau}(xR_0)')
    subplot(5,1,2)
    semilogx(Results_Matrix(6,:),Results_Matrix(1,:)*1e6,'r',tau_th,Rmax_th*1e6,'o','LineWidth',1.5);
    ylabel('R_{max} (\mum)')
    subplot(5,1,3)
    loglog(Results_Matrix(6,:),Results_Matrix(4,:)/1e6,'k','LineWidth',1.5);
    ylabel('Pressure (MPa)')
    subplot(5,1,4)
    semilogx(Results_Matrix(6,:),Results_Matrix(5,:),'c','LineWidth',1.5);
    ylabel('Velocity (m/s)')
    subplot(5,1,5)
    semilogx(Results_Matrix(6,:),Results_Matrix(2,:)*1e6,'b',tau_th,Tosc_th*1e6,'o','LineWidth',1.5);
    xlabel('Pulse duration (s)')
    ylabel('T_{osc} (\mus)')

else
    tau_th  = 0;
    Tosc_th = 0;
    Rmax_th = 0;

% Plots
    figure()
    subplot(5,1,1)
    semilogx(Results_Matrix(6,:),Results_Matrix(3,:),'g','LineWidth',1.5);
    ylabel('R|_{t=2tau} (xR0)')
    subplot(5,1,2)
    semilogx(Results_Matrix(6,:),Results_Matrix(1,:)*1e6,'r','LineWidth',1.5);
    ylabel('R_{max} (\mum)')
    subplot(5,1,3)
    loglog(Results_Matrix(6,:),Results_Matrix(4,:)/1e6,'k','LineWidth',1.5);
    ylabel('Pressure (MPa)')
    subplot(5,1,4)
    semilogx(Results_Matrix(6,:),Results_Matrix(5,:),'c','LineWidth',1.5);
    ylabel('Velocity (m/s)')
    subplot(5,1,5)
    semilogx(Results_Matrix(6,:),Results_Matrix(2,:)*1e6,'b','LineWidth',1.5);
    xlabel('Pulse duration (s)')
    ylabel('T_{osc} (\mus)')

end

end
