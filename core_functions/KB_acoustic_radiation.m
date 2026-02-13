function [tindstart,tindend,dtkb,varargout] = KB_acoustic_radiation(t,R,P,U,tind_sw,varargin)
% This function simulates acoustic radiation in the liquid around the
% bubble using the Kirkwood-Bethe hyposis See Ref:
% Liang et al. J. Fluid. Mech. 940, A5 (2022). DOI: 10.1017/jfm.2022.202
%
% Usage:
% [varargout] = KB_acoustic_radiation(t,R,P,U,tind_sw)
%
% Explicit input arguments [in SI units]:
% t     - time axis for bubble dynamics
% R     - bubble radius
% P     - pressure inside bubble
% U     - bubble wall velocity
% tind_sw - time index for selecting shock wave simulation
%
% Optional parameters
% isARexp - choice of shock wave phase. Options: 
%           1 - bubble expansion 
%           0 - bubble collapse and rebound
% winSel  - choice of selecting simulation time window. Options:
%           'cursor' - using the mouse cursor with visual interaction
%           'manual' - manually given time instants t_sw0 and t_swd
% t_sw0   - start of simulation for acoustic radiation
% t_swd   - end of simulation for acoustic radiation
% dtkb    - time step for the acoustic radiation using the KB algorithm
% p0      - ambient pressure % Update version 1.1
%
% Explicit output arguments:
% tindstart - index of time array t() where simulation of sw starts
% tindend   - index of time array t() where simulation of sw ends
%
% Optional output parameters:
% r_container - propagation distance, array size [T_KB,T_Gilmore]
% u_container - local velocity along propagating distance, array size [T_KB,T_Gilmore]
% p_container - pressure along propagating distance, array size [T_KB,T_Gilmore]
% Rtt - bubble radius interpolated on the KB time axis
% Utt - bubble wall velocity interpolated on the KB time axis
% Ptt - pressure inside the bubble interpolated on the KB time axis
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

% Notes
% Characteristics start at the bubble wall at times t(tindstart):t(tindend) with initial conditions given by the values at bubble wall
% If the choice of selecting simulation time window is set as 'cursor', then go to Figure (6) and use the mouse to click to define start time t(tindstart) and end time t(tindend)
% Remark: the R(t) curve in Fig. 6 is distorted in this graph, because it is displayed with time index. Figure 7 displays again the
% selected range of R(t) values with the correct time axis.
% 
% Some hints: for looking at sound generation at the collapse, include some time before the collapse to obtain nice p(r) curves
% For pressure wave generation during initial bubble growth, the index of time,tindstart, is automatically set to 2 and the start time is t(2)
% 
% Note that we have two time axis in simulating acoustic radiation using the KB method
% 1. tkb  = 0:dtkb:t(tindend)-t(tindstart);  % time axle for propagation of characteristics T_KB
%    dtkb - time step for the acoustic radiation using the KB algorithm - important parameter
% 2. t(tindstart:tindend);                   % time axis for bubble dynamics T_Gilmore

%% Parse optional input arguments
p = inputParser;
defaultdtkb  = 2e-12;
defaultt_sw0 = t(2);  % t(1) = 0, so we set here a non-zero minimal value
defaultt_swd = 15e-9; % 15 ns
defaultisARexp = 1;   % expansion for default
defaultwinSel  = 'cursor';
defaultp0      = 1e5;    % default value for ambient pressure 1e5 Pa

validationFcnwinSel  = @(str) any(strcmp(str,{'cursor','manuel'}));
validScalarPosNum    = @(x) isnumeric(x) && isscalar(x) && (x > 0);

addParameter(p,'isARexp',defaultisARexp);
addParameter(p,'dtkb',   defaultdtkb,validScalarPosNum);
addParameter(p,'winSel', defaultwinSel,validationFcnwinSel);
addParameter(p,'t_sw0',  defaultt_sw0,validScalarPosNum);
addParameter(p,'t_swd',  defaultt_swd,validScalarPosNum);
addParameter(p,'p0',     defaultp0);

parse(p,varargin{:});   % Parse param.-value pairs
param = p.Results;      % transfer res. to structure
clear p                 % delete the object p that is not needed

t_sw0 = param.t_sw0;
t_swd = param.t_swd;
dtkb  = param.dtkb;
isARexp=param.isARexp;
p0     = param.p0;      % ambient pressure as an input parameter

% mass density and speed of sound at ambient pressure p0
[rho, c0] = TaitEOS(p0);  
%% Constants
% p0      = 1e5;         % ambient pressure in Pa
% rho     = 998;         % density of water in ambient condiction  998kg/m3 

sigma   = 0.072583;    % surface tension 0.072583N/m
mu      = 1.046e-3;    % viscosity 0.001046 Pa.s % 7.173mPa.s for 5% dextran
B       = 3.14e8;      % parameter in Tait EOS
A       = B + p0;      % parameter in Tait EOS, A = B+p0 =3.14e8+1e5 = 3141e5;
n       = 7;           % adiabatic coeffecient in Tait EOS
fontsize = 14;         % fontsize

%% select time_index_start and time_index_end from Gilmore time axis R(t), with time window t(tindstart:tindend)

if isARexp == 1  % expansion phase
        disp('----       acoustic radiation during bubble expansion       ----')
     figure();
     plot((1:tind_sw(1)),R(1:tind_sw(1))*1e6,'r', 'LineWidth',1.5);
     xlabel('Time (index)')
     ylabel('Radius (\mum)')
     title('click to select sim. window')
     set(gca,'fontsize',fontsize)
         
         switch param.winSel
             case 'cursor'
                 disp('click on the figure to select the simulation window')
                 tclick    = ginput(2);
                 tindstart = round(tclick(1,1));
                 tindend   = round(tclick(2,1)); 
             case 'manuel'
                 % tindstart   = 1;
                 [~,tindstart] = min(abs(t-t_sw0)); %
                 [~,tindend] = min(abs(t-t_swd)); %
                 assert(t_sw0>=0);
                 % assert(t_swd<t(tind_sw(1))); 
         end

else % collapse and rebound phase
        disp('----     acoustic radiation during collapse and rebound     ----')
     figure();
     plot((tind_sw(2):tind_sw(3)), R(tind_sw(2):tind_sw(3))*1e6,'r','LineWidth',1.5);
     xlabel('Time (index)')
     ylabel('R (\mum)')
     title('click to select sim. window')
     set(gca,'fontsize',fontsize)

         switch param.winSel
             case 'cursor'
                 disp('click on the figure to select the simulation window')
                 tclick    = ginput(2);
                 tindstart = round(tclick(1,1));
                 tindend   = round(tclick(2,1));
             case 'manuel'
                 [~,tindstart] = min(abs(t-t_sw0)); %
                 [~,tindend]   = min(abs(t-t_swd));   %
                 assert(t_sw0>=t(tind_sw(2))&& t_sw0 <= t(tind_sw(3)),'error in t_sw_start');
                 assert(t_swd  > t_sw0 && t_swd <= t(tind_sw(3)),'error in t_sw_end');
         end

end

%% Plot R(t) in the selected time window t(tindstart:tindend)
figure()
plot(t(tindstart:tindend)*1e9, R(tindstart:tindend)*1e6,'LineWidth',1.5);
xlabel('Time (ns)');
ylabel('Radius (\mum)');
set(gca,'fontsize',fontsize);


%% Generate a time axis for calculating propagation of the characteristics. 
% Here we need a separate time axis, since the axis for the bubble dynamics R(t)has unequal increments.
% The time axis of each characteristic starts at zero and ends at the time given by the time window selected above.
% This second time window is set as large as t(tindstart:tindend) along time axis R(t)
% This time span can be extended, but usually much of the interesting action happens within this interval.

tkb  = 0:dtkb:t(tindend)-t(tindstart);  % time axle for propagation of characteristics [1,T_KB]

%initialize arrays for radius and pressure and velocity
r_container = zeros(length(tkb),tindend-tindstart+1); % [T_KB,T_Gilmore] array of propagation distance
p_container = zeros(size(r_container));               % [T_KB,T_Gilmore] array of pressure along propagating distance
u_container = zeros(size(r_container));               % [T_KB,T_Gilmore] array of local velocity along propagating distance

% common time axis for all characteristics
% this time axis starts again at t(tindstart), where t is the time axis 
% of the bubble motion
t1 = tkb + t(tindstart); % array t1 is used for maping the different time array for individual outgoing "characteristic" to the same array t1.

P0 = P - 2*sigma./R -4*mu*U./R;        % Pressure at bubble wall that is identical to p_bubble_wall as in gilmore_main_van_der_Waal

disp('-----                Simulating acoustic radiation         -----');
for tind = tindstart:tindend     
    
    tstart = t(tind);   %the current characteristic starts at tstart
    R0 = R(tind);       %these are the initial values for the calculation
    U0 = U(tind);
    
    % Enthalpy on the bubble wall
    H = n/(n-1)*A^(1/n)/rho*(( P0(tind) + B)^((n-1)/n)-(p0+B)^((n-1)/n));

    %Quantity y stays constant along a characteristic, which is the path of 
    %a point that moves with speed c + u outwards
    y = R0*(H + U0^2/2);
    
    % z0 contains the initial values for the differential equation
    % tt (equal to the time vector tkb)
    % and z, which contains the particle velocities and radii along the
    % characteristic, y is input parameter
    % tt = tkb  = 0:dtkb:t(tindend)-t(tindstart);  % time axle for propagation of characteristics [1,T_KB]
    z0 = [U0 R0];
    [tt,z] = ode45(@(tt,z)Kirkwood_Bethe_method(tt,z,y,'p0',p0),tkb,z0);  % tkb is the timespan [t1 tend] but it can also be a fixed time axle
    
    % we interpolate the values of r and p onto the common time axis t1,
    % because later we want to display u(t) p(r) at defined instants
    utemp = z(:,1);  % radical particle velocity along the characteristic
    rtemp = z(:,2);  % radical distance along the characteristic
    ptemp = (p0 + B)*((y./z(:,2) - z(:,1).^2/2)*(n-1)*rho/(n*(p0 + B))+1).^(n/(n-1))-B;% pressure along the characteristic
    % tstart = t(tind);   %the current characteristic starts at tstart
    % tt = tkb  = 0:dtkb:t(tindend)-t(tindstart)
    % t1 = tkb + t(tindstart);
    u_container(:,tind+1-tindstart) = interp1(tstart+tt, utemp, t1,'linear');
    r_container(:,tind+1-tindstart) = interp1(tstart+tt, rtemp, t1,'linear');
    p_container(:,tind+1-tindstart) = interp1(tstart+tt, ptemp, t1,'linear');
    
end

%For display we need the values of bubble radius and pressure interpolated
%to the new time axis t1 

Rtt = interp1(t,R,t1,'linear');
Ptt = interp1(t,P0,t1,'linear');
Utt = interp1(t,U,t1,'linear');

%We can take a look at the calculated values of r and p.
%Values for the characteristics are found in the columns.
%For display of pressure waves, we plot p versus r in a row, i.e. at the
%same time instant. This is done in the separate code kbplot.m

figure()
imagesc(t(tindstart:tindend)*1e6,t1*1e6,log10(p_container))
axis xy
colorbar
xlabel('Time along bubble dynamics (\mus)')
ylabel('Time along outgoing characteristics (\mus)')
title('pressure in the liquid, in log10 scale')
set(gca,'fontsize',12);


figure()
imagesc(t(tindstart:tindend)*1e6,t1*1e6,u_container)
axis xy
colorbar
xlabel('Time along bubble dynamics (\mus)')
ylabel('Time along outgoing characteristics (\mus)')
title('velocity in the liquid (m/s)')
set(gca,'fontsize',12);

figure()
imagesc(t(tindstart:tindend)*1e6,t1*1e6,r_container*1e6)
axis xy
colorbar
xlabel('Time along bubble dynamics (\mus)')
ylabel('Time along outgoing characteristics (\mus)')
title('progagating distance (in \mum)')
set(gca,'fontsize',12);

%% save results
varargout{1} = r_container;
varargout{2} = u_container;
varargout{3} = p_container;
varargout{4} = Rtt;
varargout{5} = Utt;
varargout{6} = Ptt;

end