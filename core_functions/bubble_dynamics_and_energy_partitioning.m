function [filename, varargout] = bubble_dynamics_and_energy_partitioning (T1,T2,R0,Rnbd,Rnc1,Rnc2,tp,varargin)
% This function use the Runge-Kutta method to numerically solve the
% extended Gilmore model and calculate the energy partitioning. See Ref:
% Liang et al. J. Fluid. Mech. 940, A5 (2022). DOI: 10.1017/jfm.2022.202

% and Section II A. in the Ref:
% Liang&Vogel,https://doi.org/10.48550/arXiv.2501.13749
%
% Usage:
% [filename, varargout] = bubble_dynamics_and_energy_partitioning (T1,T2,R0,Rnbd,Rnc1,Rnc2,tp)
%
% Explicit input arguments [in SI units]:
% T1     - first bubble oscillation time
% T2     - second bubble oscillation time
% R0     - initial bubble radius
% Rnbd   - equilibrium bubble radius after optical breakdown
% Rnc1   - equilibrium bubble radius at the first bubble collapse
% Rnc2   - equilibrium bubble radius at the second bubble collapse
% tp     - pulse duration
%
% Optional parameters
% GS     - choice of the Gilmore Solver (GS). Options:
%           'adaptive' - fast, ideal for bubble dynamics and energy partitioning
%           'maxstepsize'-adaptive method with maximum stepsize control,
%           ideal for simulating acoustic emission that requires finer
%           temporal grids, especially for emission at bubble collapse
%           'fixedstepsize'-fixed stepsize
% MaxStepSize - tuning max. step size when GilmoreSolver = 'maxstepsize'is selected
% Stepsize    - tuning step size when GilmoreSolver = 'fixedstepsize'is selected
% EP     - choice of calculating Energy Partitioning (EP). Options:
%                1  - performing calculation of energy partitioning
%                0  - not performed
% tend   - time instant that simulation ends
%
% Optional parameters used in the function extended_gilmore
% RNP    - radius of extra hard core, e.g. nanoparticles
% tRmax3 - time instant at Rmax3
% Twall  - temperature at the bubble wall
% JSC    - choice of the jump-start condition (JSC). Options:
%                'no' - no jump-start condition
%                'first-order'- first-order approximation for inertial
%                confined condition
%                'second-order'- second-order approximation for inertial
%                confined condition
%                'general' - second-order approximation for general
%                condition (default)
%
% Output arguments:
% filename  - name of the excel file to store the results
% varargout - varable-length argument output
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
defaultJSC    = 'general'; % default value for jump-start-condition (JSC)
defaultRNP    = 0;         % default value for the radius of hard core, e.g. nanoparticles
defaulttRmax3 = 1000;      % default value for tRmax, set a large value as "infinite", in seconds
defaultTwall  = 293;       % default value for the temperature at bubble wall,in Kelvin
defaultRnc3   = Rnc2;      % default value for Rncoll3
defaultisIsotherm  = 1;    % default value for the kappa value for late oscillations t > tRmax3, 1 for isothermal, 0 for adiabatic

defaultGS     = 'adaptive';   % default value for Gilmore Solver (GS)
defaultMaxStepSize = 200e-12; % default value for MaxStepSize; MaxStepSize = 100e-12; % 200ps for Rmax = 5.4 um, 1ns for Rmax >= 100 um and 100 ps for Rmax < 1 um
defaultStepSize    = 20e-12;  % default value for StepSize 20 ps
defaultEP     = 1;            % default value for EnergyPartitioning (EP)
defaulttend   = (T1+T2)*1.5;  % default value for tend
defaultp0     = 1e5;          % default value for ambient pressure 1e5 Pa

validationFcnJSC  = @(str) any(strcmp(str,{'no','first-order','second-order','general'}));
validationFcnGS   = @(str) any(strcmp(str,{'adaptive','maxstepsize','fixedstepsize'}));
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);     % positive scalar number
validScalarnonNegNum = @(x) isnumeric(x) && isscalar(x) && (x >= 0); % non-negative scalar number

addParameter(p,'JSC',defaultJSC,validationFcnJSC);
addParameter(p,'RNP',defaultRNP,validScalarnonNegNum);
addParameter(p,'tRmax3',defaulttRmax3,validScalarPosNum);
addParameter(p,'Twall', defaultTwall,validScalarPosNum);
addParameter(p,'Rnc3', defaultRnc3);
addParameter(p,'isIsotherm', defaultisIsotherm);

addParameter(p,'GS',defaultGS,validationFcnGS);
addParameter(p,'MaxStepSize', defaultMaxStepSize,validScalarPosNum);
addParameter(p,'StepSize', defaultStepSize,validScalarPosNum);
addParameter(p,'EP', defaultEP,validScalarnonNegNum);
addParameter(p,'tend',defaulttend,validScalarPosNum);
addParameter(p,'p0',defaultp0);


parse(p,varargin{:});   % Parse param.-value pairs
param = p.Results;      % transfer res. to structure
clear p                 % delete the object p that is not needed

RNP = param.RNP;
p0     = param.p0;      % ambient pressure as an input parameter

% mass density and speed of sound at ambient pressure p0
[rho, c0] = TaitEOS(p0);  

%% Constants
% p0      = 1e5;         % ambient pressure in Pa
% rho     = 998;         % density of water in ambient condiction  998kg/m3
% c0      = 1483;        % speed of sound 1483m/s
pv      = 2330;        % vapor pressure at 20 degree in Pa
sigma   = 0.072583;    % surface tension 0.072583N/m
kappa   = 4/3;         % adiabatic exponent 4/3 for H2O vapor; 5/3 for monoatomic gas; 7/5 for diatomic gas 
mu      = 1.046e-3;    % viscosity 0.001046 Pa.s % 7.173mPa.s for 5% dextran
B       = 3.14e8;      % parameter in Tait EOS
A       = B+p0;        % parameter in Tait EOS, A = B+p0 =3.14e8+1e5 = 3141e5;
n       = 7;           % adiabatic coeffecient in Tait EOS
Cp      = 4187;        % J/(K kg) at 20 °C
L_w     = 2256*1e3;    % J/kg latent heat at  100 °C
rho_vap = 0.761;       % kg/m^3 mass dnesity of vapor at 20 degree celsius and 1 bar
fonts   = 14;          % font size for graphics

%% Input parameters -----
%switch Rn(t) at maximum bubble expansions
tRmax  = T1/2;        % at tRmax1 Rnt is switched from Rnbd to Rnc1
tRmax2 = T1+T2/2;     % at tRmax2 Rnt is switched from Rnc1 to Rnc2

% start and end of calculation, time span for simulation is set as [t1 t2]
t1 = 0;               % starting time for the calculation
t2 = param.tend;      % end time for the calculation, default value (T1+T2)*1.5

% check input times
assert(tp>0&&tp<1,'Error in Pulse duration Input');
assert(T1>0&&T1<1,'Error in Tosc1 Input');
assert(T2>0&&T2<1,'Error in Tosc2 Input');

%% name of the output excel file in full path
toolboxPath = fileparts(fileparts(mfilename('fullpath'))); % Get the toolbox path
filename = ['Report_T1_' num2str(T1/1e-6) 'us_T2_' num2str(T2/1e-6) 'us_pulse_' num2str(tp/1e-9) 'ns.xlsx'];
filename = fullfile(toolboxPath,'results',filename);

%% initial and equilibrium radii of the bubble
Rn = Rnbd ;   % equilibrium radius for breakdown, Rnbd
Rn2= Rnc1;    % equilibrium radius for 1st collapse when t>t_Rmax  Rncoll
Rn3= Rnc2;    % equilibrium radius for 2nd collapse when t>t_Rmax2 
           
%% initial condition and run ODE Runge-Kutta 

x0 = [R0 0];  % initial value of  R and U, here R = R0 and U =0;

% solving Gilmore model using the adaptive 45 Runge-Kutta method
switch param.GS        % choice of the GilmoreSolver
    case 'adaptive'    % Adaptive 45 Runge-Kutta method
    disp('Adaptive 45 Runge-Kutta method');  
    opts = odeset('RelTol',1e-8,'AbsTol',1e-9,'NormControl','off'); % Default: RelTol = 1e-3; AbsTol = 1e-6;
    [t,x] = ode45(@(t,x)extended_gilmore(t,x,R0,Rn,Rn2,Rn3,tRmax,tRmax2,tp,'JSC',param.JSC,'RNP',RNP,'tRmax3',param.tRmax3,'Twall',param.Twall,'isIsotherm',param.isIsotherm,'Rnc3',param.Rnc3,'p0',p0),[t1 t2],x0,opts);
    
    case 'maxstepsize' % Adaptive 45 Runge-Kutta method with max. stepsize ideal for simulating shock wave emission
    disp(['Adaptive 45 Runge-Kutta method with max. stepsize control of ' num2str(param.MaxStepSize)]);
    opts = odeset('RelTol',1e-5,'AbsTol',1e-8,'MaxStep',param.MaxStepSize);
    [t,x] = ode45(@(t,x)extended_gilmore(t,x,R0,Rn,Rn2,Rn3,tRmax,tRmax2,tp,'JSC',param.JSC,'RNP',RNP,'tRmax3',param.tRmax3,'Twall',param.Twall,'isIsotherm',param.isIsotherm,'Rnc3',param.Rnc3,'p0',p0),[t1 t2],x0,opts);
    
    case 'fixedstepsize'  % Runge-Kutta method with fixed stepsize, slow
    disp(['Adaptive 45 Runge-Kutta method with fixed stepsize control of' num2str(param.StepSize)]);
    timespan = t1:param.StepSize:t2;
    [t,x] = ode45(@(t,x)extended_gilmore(t,x,R0,Rn,Rn2,Rn3,tRmax,tRmax2,tp,'JSC',param.JSC,'RNP',RNP,'tRmax3',param.tRmax3,'Twall',param.Twall,'isIsotherm',param.isIsotherm,'Rnc3',param.Rnc3,'p0',p0),timespan,x0);   
end

%write R and U in separate vectors
R = x(:,1);  % bubble radius
U = x(:,2);  % bubble wall speed


%% Plots
%plot the radius-time curve
gcf = figure();
set(gcf,'units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
plot(t*1e6,R*1e6,'LineWidth',1.5);
xlim([-0.05 1.00]*max(t)*1e6);
xlabel('Time (\mus)')
ylabel('Bubble radius (\mum)')
set(gca,'fontsize',fonts)
%% Calculate the maximum radius Rmax
[Rmax,mt] = max(R);
index_Rmax= mt;
t_Rmax    = t(mt);
Rratio    = Rmax/R0;

%% Calculate the Rmin and the 1st osc. time Tosc1

[Rmin,nt]=min(R(mt:length(R))); % find the Rmin
index_Rmin=nt+mt-1;             % index of Rmin.

[pks,locs] = findpeaks(R(1:index_Rmin),t(1:index_Rmin));
if length(pks)==1 % to check if the Rmin. is the first collapse
    tRmin = t(index_Rmin); 
else
    [pks2,locs2] = findpeaks(R(1:index_Rmin));
    [Rmin,index_Rmin] = min(R(locs2(1):locs2(2)));% find the minima inbetween the 1st peak and the 2nd peak
    index_Rmin = index_Rmin+locs2(1)-1;
    tRmin = t(index_Rmin); 
end

% Tosc_2 and Rmax_2
[Rmax2,index_Rmax2] = max(R(index_Rmin:end));
index_Rmax2 = index_Rmax2+index_Rmin-1;
[Rmin2,index_Rmin2] = min(R(index_Rmax2:end)); % from t_Rmax2 to the end, search for Tosc2
index_Rmin2 = index_Rmax2+index_Rmin2-1;

tau   = t(index_Rmin)*1e9;  % osc. time in ns
Tosc1 = tau;                % 1st osc. time in ns
Tosc2 = t(index_Rmin2)*1e9; % 2nd osc. time in ns

time_expansion = t_Rmax*1e9;% in ns expansion time
time_collapse  = ((t(index_Rmin)-t_Rmax))*1e9; %in ns collapse time
Tcoll_Texp_ratio = (t(index_Rmin)-t_Rmax)/t_Rmax;
title(['\tau_o_s_c = ' num2str(tau,3) ' ns, R_m_a_x = ' num2str(Rmax*1e6,3) ' um, R_m_i_n = ' num2str(Rmin*1e9,3) ' nm, R_m_a_x_2 = ' num2str(Rmax2*1e6,3) 'um'],'fontsize',12)

%plot the velocity-time curve
gcf = figure();
set(gcf,'units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
plot(t*1e6,U,'LineWidth',1.5);
xlim([-0.05 1.00]*max(t)*1e6);
set(gca,'fontsize',fonts)
xlabel('Time (\mus)')
ylabel('Velocity (m/s)')
title('Bubble wall velocity (m/s) ');


%% Latent heat of vaporization
Rvap_bd    = R0*(rho/rho_vap)^(1/3); % vapor radius after bd at 20 degree celsius and 1 bar
Rvap_Rmax1 = (pv/p0)^(1/3)*Rmax;     % vapor radius at Rmax1 at 20 degree celsius and 1 bar

if Rn2 > Rn3 % Info on Tosc1, Tosc2, Tosc3 are available, e.g. JFM paper
    Rvap_coll  = (Rn2^3-Rn3^3)^(1/3);% vapor radius at Rmin1 at 20 degree celsius and 1 bar
else         % Info on Tosc1, Tosc2 are available but on Tosc3 is not available, e.g. 560 ps, Tosc1, Tosc2 (Rmax) data
    Rvap_coll = Rn2;                 % vapor radius at Rmin1 at 20 degree celsius and 1 bar
end

Rvap_Rmax2 = (pv/p0)^(1/3)*Rmax2;    % vapor radius at Rmax2 at 20 degree celsius and 1 bar

E_v_tot    = 4/3*pi*R0^3*rho*(Cp*(100-20)+L_w); % total vapor energy in J
E_v_tot_nJ = E_v_tot*1e9;                       % total vapor energy in nJ
E_v_Rmax1  = Rvap_Rmax1^3/Rvap_bd^3*E_v_tot;    % vapor energy in J at Rmax1
E_v_coll   = Rvap_coll^3/Rvap_bd^3*E_v_tot;     % vapor energy in J at collapse
E_v_Rmax2  = Rvap_Rmax2^3/Rvap_bd^3*E_v_tot;    % vapor energy in J at Rmax2

Delta_E_v_exp = E_v_tot - E_v_Rmax1;            % dissipation of vapor during expasion in J
Delta_E_v_coll_reb = E_v_Rmax1 - E_v_Rmax2;     % dissipation of vapor during collapse and rebound in J

E_v_t       = E_v_tot/tp*(t/2-tp/(2*pi)*sin(pi*t/tp));
E_v_t       = E_v_t.*(t<=2*tp) + E_v_tot*(t>2*tp & t<=tRmax) + ...
              E_v_Rmax1*(t>tRmax & t<=tRmin) + E_v_coll*(t>tRmin & t<=tRmax2) + E_v_Rmax2*(t>tRmax2);
Delta_E_v_t = 0*(t<=tRmax) + (E_v_tot-E_v_Rmax1)*(t>tRmax & t<=tRmin) +...
             (E_v_tot-E_v_coll)*(t>tRmin & t<=tRmax2) + (E_v_tot-E_v_Rmax2)*(t>tRmax2);

%% the equilibrium bubble radius Rn(t)
Rnt = (R0^3 + (Rn^3-R0^3)/2/tp*(t-tp/pi*sin(pi/tp*t))).^(1/3);
Rnt = Rnt.*(t<=2*tp) + Rn*(t>2*tp & t<=tRmax) + Rn2*(t>tRmax & t<=tRmax2)+Rn3*(t>tRmax2);   %

%% van der Waals hard core Rvdw considering extra hard core as NPs
% Here van der waals' hard core is implemented
if RNP > 0
    Rvdw = RNP*(t<tRmax)+(((1/9)^3*(Rnt.^3-RNP^3)+RNP^3).^(1/3)).*(t>=tRmax); % before Rvdw = RNP
    assert(R0 > RNP);  % radius of the initial hot vapor size should be larger than the nanoparticle size.
else
    Rvdw=0*(t<tRmax)+Rnt/9.*(t>=tRmax);
end

gcf = figure();
set(gcf,'units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
plot(t*1e6,Rnt*1e6,t*1e6,Rvdw*1e6, t*1e6,R*1e6,'LineWidth',1.5);
xlim([-0.05 1.00]*max(t)*1e6);
ylim([-0.05 1.1]*max(R)*1e6);
set(gca,'fontsize',fonts)
xlabel('Time (\mus)')
ylabel('Radius (\mum)')
legend('R_n(t)','R_v_d_w(t)','R(t)');
title('R(t), R_n(t) and R_v_d_w(t)');


pn = p0 + 2*sigma./Rnt;                              % gas pressure at the current value of the equilibrium radius, Rnt
pg = pn.*((Rnt.^3-Rvdw.^3)./(R.^3-Rvdw.^3)).^kappa;  % gas pressure in the bubble at the current radius R = R(t) with van der Waals


pgas_exp_max = max(pg(1:index_Rmax));           % in unit Pa
U_exp_max    = max(U(1:index_Rmax));            % in m/s
pgas_coll_max= max(pg(index_Rmax:index_Rmax2)); % in unit Pa
U_coll_min   = min(U(index_Rmax:index_Rmax2));  % in m/s
U_reb_max    = max(U(index_Rmax:index_Rmax2));  % in m/s

pgas_max1_up   = (p0+2*sigma/Rn)*(Rn/Rmax)^(3*kappa);                           % Rvdw =0 in Pa
pgas_max1_down = (p0+2*sigma/Rn2)*((Rn2^3-(Rn2/9)^3)/(Rmax^3-(Rn2/9)^3))^kappa; % Rvdw = Rn/9 in Pa

p_bubble_wall = pg - 2*sigma./R -4*mu*U./R;            % pressure on the bubble wall, which meets the Kickwood-Bethe hypothesis and propagate outward
H_enthaphy    = n/(n-1)*A^(1/n)/rho*((p_bubble_wall+B).^((n-1)/n)-(p0+B)^((n-1)/n)); % Enthaphy in the liquid
C_bubble_wall = (c0^2+(n-1)*H_enthaphy).^0.5;          % speed of sound at the bubble wall
norm_density  = ((p_bubble_wall+B)/(p0+B)).^(1/n);

% find the region where SW is able to be produced (identify bubble wall pressure > p0 = 1e5 Pa)
tind_sw1 = find(p_bubble_wall(1:index_Rmin)<p0,1,'first')-1;
tind_sw2 = find(p_bubble_wall(1:index_Rmin)<p0,1,'last')+1;
tind_sw3 = find(p_bubble_wall(index_Rmin:end)<p0,1,'first')+index_Rmin-2;

% key time index that defines the windown for shock wave simulation
tind_sw = [tind_sw1 tind_sw2 tind_sw3];

%plot the pressure inside the bubble in a semilog scale
gcf = figure();
set(gcf,'units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
semilogy(t*1e6,pg/1e6,'LineWidth',1.5);
xlim([-0.05 1.00]*max(t)*1e6);
set(gca,'fontsize',fonts);
xlabel('Time (\mus)')
ylabel('Pressure (MPa)')
title('Pressure inside the bubble');
hold on;
semilogy([t(tind_sw1) t(tind_sw2) t(tind_sw3)]*1e6,[pg(tind_sw1) pg(tind_sw2) pg(tind_sw3)]/1e6,'ro');
hold off;

p_surf = 2*sigma./R;
p_visc = 4*mu./R.*U;

%plot all pressures
gcf=figure();
set(gcf,'units','normalized','outerposition',[0.05 0.05 0.9 0.9]);

subplot(2,2,1);
semilogy(t*1e6,pg/1e6,'LineWidth',1.5);
xlim([-0.02 1.00]*max(t)*1e6);
set(gca,'fontsize',fonts);
xlabel('Time (\mus)');
ylabel('Pressure (MPa)');
legend('p_g_a_s');
txt1 = ['\leftarrow ' num2str(max(pg(1:index_Rmax)/1e6),'%.1f') ' MPa'];
text(0,max(pg(1:index_Rmax)/1e6),txt1,'HorizontalAlignment','left');

subplot(2,2,2);
plot(t*1e6,(p_bubble_wall)/1e6,'LineWidth',1.5);
xlim([-0.02 1.00]*max(t)*1e6);
ylim([-0.05 1.1]*max((p_bubble_wall)/1e6));
set(gca,'fontsize',fonts);
xlabel('Time (\mus)');
ylabel('Pressure (MPa)');
legend('p_w_a_l_l');
txt1 = ['\leftarrow ' num2str(max(p_bubble_wall(1:index_Rmax)/1e6),'%.1f') ' MPa'];
text(0,max(p_bubble_wall(1:index_Rmax)/1e6),txt1,'HorizontalAlignment','left');

subplot(2,2,3);
semilogy(t*1e6,p_surf/1e6,'LineWidth',1.5);
xlim([-0.02 1.00]*max(t)*1e6);
set(gca,'fontsize',fonts);
xlabel('Time (\mus)');
ylabel('Pressure (MPa)');
legend('p_s_u_r_f');
txt1 = ['\leftarrow ' num2str(max(p_surf(1:index_Rmax)/1e6),'%.1f') ' MPa'];
text(0,max(p_surf(1:index_Rmax)/1e6),txt1,'HorizontalAlignment','left');

subplot(2,2,4);
plot(t*1e6,p_visc/1e6,'LineWidth',1.5);
xlim([-0.02 1.00]*max(t)*1e6);
set(gca,'fontsize',fonts);
xlabel('Time (\mus)');
ylabel('Pressure (MPa)');
legend('p_v_i_s_c');
txt1 = ['\leftarrow ' num2str(max(p_visc(1:index_Rmax)/1e6),'%.1f') ' MPa'];
text(0,max(p_visc(1:index_Rmax)/1e6),txt1,'HorizontalAlignment','left');

%% Energy partitioning
if param.EP == 1
    disp('---              Energy partitioning is switched on           ---');
% Internal energy U_gas
U_gas = 4*pi/(3*(kappa-1))*(pg.*(R.^3-Rvdw.^3)-pg(1)*R(1)^3); %Note: van der Waal's hard core should be taken into account

% Condensation of internal energy, Delta_U_int,cond
Delta_U_int_cond_Rmax1 = 4*pi/(3*(kappa-1))*((p0+2*sigma/Rn)*(Rn/Rmax)^(3*kappa)*Rmax^3-(p0+2*sigma/Rn2)*(Rn2/Rmax)^(3*kappa)*Rmax^3); % Rn = Rnbd, Rn2 = Rncoll
Delta_U_int_cond_Rmax2 = 4*pi/(3*(kappa-1))*Rmax2^3*((p0+2*sigma/Rn2)*(Rn2/Rmax2)^(3*kappa)-(p0+2*sigma/Rn3)*(Rn3/Rmax2)^(3*kappa));   % Rn = Rnbd, Rn2 = Rncoll,Rn3=Rncoll2
Delta_U_int_cond       = 0*(t<=tRmax) + (Delta_U_int_cond_Rmax1)*(t>tRmax & t<=tRmax2) + (Delta_U_int_cond_Rmax1+Delta_U_int_cond_Rmax2)*(t>tRmax2);

% Work done by the gas, numerical integral(4*pi*R^2*U*pgas)dt
Work_gas_integrand  = 4*pi*((R.*pg).*R).*U;  % 4*pi*R.^2.*U.*pg
Work_visc_integrand = 16*pi*mu*R.*U.^2;
Work_surf_integrand = 8*pi*sigma*R.*U;
Work_pot_integrand  = 8*pi*sigma*R.*U+4*pi*p0*R.*U.*R;

Work_gas_integral   = ones(length(t),1);
Work_visc_integral  = ones(length(t),1);
Work_surf_integral  = ones(length(t),1);
Work_pot_integral   = ones(length(t),1);

Work_gas_integral(1)  = 0;
Work_visc_integral(1) = 0;
Work_surf_integral(1) = 0;
Work_pot_integral(1)  = 0;

for i=2:length(t)
    Work_gas_integral(i)  = trapz(t(1:i),Work_gas_integrand(1:i));
    Work_visc_integral(i) = trapz(t(1:i),Work_visc_integrand(1:i));
    Work_surf_integral(i) = trapz(t(1:i),Work_surf_integrand(1:i));
    Work_pot_integral(i)  = trapz(t(1:i),Work_pot_integrand(1:i));
end

E_total_num = U_gas + Work_gas_integral + ...
              E_v_t + Delta_E_v_t + Delta_U_int_cond;

gcf = figure();
set(gcf,'units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
plot(t*1e6,E_total_num*1e9,t*1e6,U_gas*1e9,t*1e6,Delta_U_int_cond*1e9,t*1e6,Work_gas_integral*1e9,t*1e6,E_v_t*1e9,t*1e6,Delta_E_v_t*1e9, t*1e6,Work_visc_integral*1e9, t*1e6,Work_surf_integral*1e9,t*1e6,Work_pot_integral*1e9,'LineWidth',1.5);
xlim([-0.05 1.0]*max(t)*1e6);
ylim([-0.05 1.05]*max(E_total_num)*1e9);
set(gca,'fontsize',fonts);
xlabel('Time (\mus)');
ylabel('Energy (nJ)');
legend('E_t_o_t_a_l','U_{int}','\DeltaU_{int,cond}','W_g_a_s','E_{vap}','\DeltaE_{vap}','W_v_i_s_c','W_s_u_r_f','E_p_o_t');
title('energy partitioning');

%% Total energy
% total energy - Etotal
[max_E_total_num,index_max_E_total_num] = max(E_total_num);
max_E_total_num_nJ = max_E_total_num*1e9;% in nJ
time_E_total_num_us = t(index_max_E_total_num)*1e6;%in us
txt1 = ['\downarrow E_t_o_t_a_l = ' num2str(max_E_total_num_nJ*1e3,'%.1f') ' pJ'];
text(time_E_total_num_us,max_E_total_num_nJ,txt1,'VerticalAlignment','bottom');

%% Energy partitioning at breakdown 2023.01.20
E_v_bd_pJ   = E_v_tot_nJ*1e3;                   % vaporization energy after breakdown in pJ
U_int_bd_pJ = max_E_total_num_nJ*1e3-E_v_bd_pJ; % internal energy after breakdown in pJ
E_v_bd_per  = E_v_bd_pJ/(max_E_total_num_nJ*1e3)*100;    % vaporization energy after breakdown in %
U_int_bd_per = U_int_bd_pJ/(max_E_total_num_nJ*1e3)*100; % internal energy after breakdown in %

%% Energy partitioning at tmax1
% max. Wgas - Wgas,max
[max_Work_gas_integral,index_max_Work_gas_integral] = max(Work_gas_integral);
max_Work_gas_integral_nJ = max_Work_gas_integral*1e9;% in nJ
time_max_Work_gas_integral_us = t(index_max_Work_gas_integral)*1e6;%in us
txt2 = ['\uparrow W_g_a_s = ' num2str(max_Work_gas_integral_nJ*1e3,'%.1f') ' pJ'];
text(time_max_Work_gas_integral_us,max_Work_gas_integral_nJ,txt2,'VerticalAlignment','top');

% max. Epotential 
[max_Work_pot_integral,index_Work_pot_integral] = max(Work_pot_integral);
max_Work_pot_integral_nJ = max_Work_pot_integral*1e9;% in nJ
index_Work_pot_integral_us = t(index_Work_pot_integral)*1e6;%in us
txt1 = ['\downarrow E_p_o_t = ' num2str(max_Work_pot_integral_nJ*1e3,'%.1f') ' pJ'];
text(index_Work_pot_integral_us,max_Work_pot_integral_nJ,txt1,'VerticalAlignment','bottom');

% max. Surf.
[max_Work_surf_integral,index_Work_surf_integral] = max(Work_surf_integral);
max_Work_surf_integral_nJ = max_Work_surf_integral*1e9;% in nJ
index_Work_surf_integral_us = t(index_Work_surf_integral)*1e6;%in us
txt1 = ['\downarrow W_s_u_r_f = ' num2str(max_Work_surf_integral_nJ*1e3,'%.1f') ' pJ'];
text(index_Work_surf_integral_us,max_Work_surf_integral_nJ,txt1,'VerticalAlignment','bottom');

% max. W_stat.
max_Work_stat_nJ = max_Work_pot_integral_nJ-max_Work_surf_integral_nJ;
txt1 = ['W_s_t_a_t = ' num2str((max_Work_stat_nJ)*1e3,'%.1f') ' pJ'];
text(index_Work_surf_integral_us,(max_Work_surf_integral_nJ+max_Work_pot_integral_nJ)/2,txt1,'VerticalAlignment','middle');

% max. Visc.
max_Work_visc_integral = Work_visc_integral(index_Work_surf_integral);
max_Work_visc_integral_nJ = max_Work_visc_integral*1e9;% in nJ
time_max_Work_visc_integral_us = t(index_Work_surf_integral)*1e6;%in us
txt2 = ['\uparrow W_v_i_s_c = ' num2str(max_Work_visc_integral_nJ*1e3,'%.1f') ' pJ'];
text(time_max_Work_visc_integral_us,max_Work_visc_integral_nJ,txt2,'VerticalAlignment','top');

% max. Shock Wave
E_SW1_nJ = (max_Work_gas_integral_nJ-max_Work_pot_integral_nJ-max_Work_visc_integral_nJ); % in nJ

txt1 = ['E_S_W = ' num2str(E_SW1_nJ*1e3,'%.1f') ' pJ'];
text(index_Work_surf_integral_us,(max_Work_gas_integral_nJ+max_Work_pot_integral_nJ)/2,txt1,'VerticalAlignment','middle');

% internal energy at tmax1, Ugas,max1 = pv*Rmax^3 - p0R0^3
U_gas_max1 = 4*pi/(3*(kappa-1))*(pv*Rmax^3-pg(1)*R0^3);
U_gas_max1_nJ = U_gas_max1*1e9; % in nJ
txt1 = ['U_g_a_s_,_m_a_x = ' num2str((U_gas_max1_nJ)*1e3,'%.1f') ' pJ \rightarrow'];
text(index_Work_surf_integral_us,U_gas_max1_nJ,txt1,'HorizontalAlignment','right');

% internal energy at t>=tmax1, Ugas,coll
U_gas_coll = 4*pi/(3*(kappa-1))*((p0+2*sigma/Rn2)*(Rn2/Rmax)^(3*kappa)*Rmax^3-pg(1)*R0^3); % Rn2 = Rncoll
U_gas_coll_nJ = U_gas_coll*1e9; % in nJ
txt1 = ['\uparrow U_g_a_s_,_c_o_l_l = ' num2str((U_gas_coll_nJ)*1e3,'%.1f') ' pJ'];
text(index_Work_surf_integral_us,U_gas_coll_nJ,txt1,'VerticalAlignment','top');

% internal energy at t<=tmax, Ugas,max1
U_gas_max1_up    = 4*pi/(3*(kappa-1))*((p0+2*sigma/Rn)*(Rn/Rmax)^(3*kappa)*Rmax^3-pg(1)*R0^3); % Rn2 = Rncoll1
U_gas_max1_up_nJ = U_gas_max1_up*1e9; % in nJ

% Vapor energy at Rmax
E_v_Rmax1_nJ     = E_v_Rmax1*1e9; % in nJ
% Vapor energy at Rmin
E_v_coll_nJ      = E_v_coll*1e9;  % in nJ
% Vapor energy at Rmax2
E_v_Rmax2_nJ     = E_v_Rmax2*1e9; % in nJ

% Condensation of vapor energy at expansion
Delta_E_v_exp_nJ = (E_v_tot-E_v_Rmax1)*1e9;   % in nJ
% Condensation of vapor energy at collapse
Delta_E_v_coll_nJ = (E_v_Rmax1-E_v_coll)*1e9; % in nJ
% Condensation of vapor energy at rebound
Delta_E_v_reb_nJ = (E_v_coll-E_v_Rmax2)*1e9;  % in nJ
% Condensation of vapor energy for residual
Delta_E_v_resid_nJ = (E_v_Rmax2)*1e9;         % in nJ

% bubble energy at Rmax = Epot + Ugas + Ev (bubble energy = potential energy + internal energy + vapor energy)
E_B_max1_nJ = U_gas_max1_nJ + max_Work_pot_integral_nJ + E_v_Rmax1_nJ;  % 

% Condensation dissipation of internal energy at Rmax1 (the 1st osc).

Delta_U_int_cond_Rmax1_nJ = Delta_U_int_cond_Rmax1 * 1e9; % in nJ
txt1 = ['\leftarrow E_c_o_n_d = ' num2str((Delta_U_int_cond_Rmax1_nJ)*1e3,'%.1f') ' pJ'];
text(index_Work_surf_integral_us,(max_E_total_num_nJ+max_Work_gas_integral_nJ)/2,txt1,'HorizontalAlignment','left');

% Condensation dissipation of internal energy at expansion
Delta_U_int_cond_exp_nJ = Delta_U_int_cond_Rmax1_nJ + U_gas_coll_nJ - U_gas_max1_nJ;% in nJ

% Condensation dissipation of internal energy and vapor energy at expansion
E_cond_exp_nJ = Delta_E_v_exp_nJ + Delta_U_int_cond_exp_nJ;

% Condensation dissipation of internal energy at collapse
Delta_U_int_cond_coll_nJ = U_gas_max1_nJ - U_gas_coll_nJ;

% Condensation dissipation of internal energy and vapor energy at collapse
E_cond_coll_nJ = Delta_E_v_coll_nJ + Delta_U_int_cond_coll_nJ;

%% Energy partitioning at Rmin1
% intermal energy at Rmin1
U_gas_min = U_gas(index_Rmin);
U_gas_min_nJ = U_gas_min*1e9; % in nJ
time_Rmin_us = t(index_Rmin)*1e6; % in us
txt1 = ['\leftarrow U_g_a_s_,_m_i_n = ' num2str((U_gas_min_nJ)*1e3,'%.1f') ' pJ'];
text(time_Rmin_us,U_gas_min_nJ,txt1,'HorizontalAlignment','left');

% viscosity at Rmin1
W_visc_min = Work_visc_integral(index_Rmin)-max_Work_visc_integral;
W_visc_min_nJ = W_visc_min*1e9; % in nJ
txt1 = ['\leftarrow W_v_i_s_c_,_m_i_n = ' num2str((W_visc_min_nJ)*1e3,'%.1f') ' pJ'];
text(time_Rmin_us,Work_visc_integral(index_Rmin)*1e9,txt1,'HorizontalAlignment','left');

% compression in liquid at Rmin1
E_comp_coll = max_Work_pot_integral+U_gas_coll-W_visc_min-U_gas_min;
E_comp_coll_nJ = E_comp_coll*1e9;  % in nJ
txt1 = ['E_c_o_m_p_,_c_o_l_l = ' num2str((E_comp_coll_nJ)*1e3,'%.1f') ' pJ'];
text(time_Rmin_us,U_gas_min_nJ+Work_visc_integral(index_Rmin)*1e9,txt1,'HorizontalAlignment','left');

%% Energy partitioning at Rmax2
% viscosity at Rmax2
W_visc_max2 = Work_visc_integral(index_Rmax2)-Work_visc_integral(index_Rmin);
W_visc_max2_nJ = W_visc_max2*1e9; % in nJ
time_Rmax2_us = t(index_Rmax2)*1e6; % in us
txt1 = ['\downarrow W_v_i_s_c_,_m_a_x_2 = ' num2str((W_visc_max2_nJ)*1e3,'%.2f') ' pJ'];
text(time_Rmax2_us,Work_visc_integral(index_Rmax2)*1e9,txt1,'VerticalAlignment','bottom');

% potential at Rmax2
E_pot_max2 = Work_pot_integral(index_Rmax2);
E_pot_max2_nJ = E_pot_max2*1e9; % in nJ
txt1 = ['\downarrow E_p_o_t_,_m_a_x_2 = ' num2str((E_pot_max2_nJ)*1e3,'%.2f') ' pJ'];
text(time_Rmax2_us,E_pot_max2_nJ,txt1,'VerticalAlignment','bottom');

% W_Surf. at Rmax2
W_surf_max2_nJ= Work_surf_integral(index_Rmax2)*1e9; % in nJ

% W_stat. at Rmax2
W_stat_max2_nJ = E_pot_max2_nJ-W_surf_max2_nJ; % in nJ

% condensation dissipation in the 2nd Osc.
E_cond_2osc_nJ = Delta_U_int_cond_Rmax2 * 1e9; % in nJ

% internal energy at tmax2, Ugas,max2 = 4pi*?pv*Rmax2^3 - p0R0^3?
U_gas_max2 = 4*pi/(3*(kappa-1))*(pv*Rmax2^3-pg(1)*R0^3);    % internal energy at Rmax2
U_gas_max2_nJ = U_gas_max2*1e9;                             % in nJ

if E_cond_2osc_nJ==0
 U_gas_max2    = U_gas(index_Rmax2);
 U_gas_max2_nJ = U_gas_max2*1e9; % if no condensation is considered in Rmax2, 
                                         % the internal energy of gas should not be calculated by pv = 2330 Pa but read from Ugas(t) curve
                                         % U_gas_max2_nJ updated
end
txt1 = ['\uparrow U_g_a_s_,_m_a_x_2 = ' num2str((U_gas_max2_nJ)*1e3,'%.2f') ' pJ'];
text(time_Rmax2_us,U_gas_max2_nJ,txt1,'VerticalAlignment','top');

% bubble energy at Rmax2
E_B_max2_nJ =W_stat_max2_nJ + W_surf_max2_nJ+U_gas_max2_nJ+E_v_Rmax2_nJ;

% internal energy at t<=tmax2, Ugas,max2
U_gas_max2_up    = 4*pi/(3*(kappa-1))*((p0+2*sigma/Rn2)*(Rn2/Rmax2)^(3*kappa)*Rmax2^3-pg(1)*R0^3); % Rn2 = Rncoll1
U_gas_max2_up_nJ = U_gas_max2_up*1e9; % in nJ

% internal energy at t>=tmax2, Ugas,max2
U_gas_max2_down = 4*pi/(3*(kappa-1))*((p0+2*sigma/Rn3)*(Rn3/Rmax2)^(3*kappa)*Rmax2^3-pg(1)*R0^3); % Rn3 = Rncoll2
U_gas_max2_down_nJ = U_gas_max2_down*1e9; % in nJ

Delta_U_int_cond_reb_nJ = U_gas_max2_up_nJ - U_gas_max2_nJ;    % condensation dissipation in rebound   (tmin1<=t<=tmax2)
Delta_U_int_cond_resid_nJ = U_gas_max2_nJ - U_gas_max2_down_nJ; % condensation dissipation in afterbounces (t>tmax2)

if E_cond_2osc_nJ == 0  % if Rn2 = Rn3, no condensation is considered!
    Delta_U_int_cond_reb_nJ = 0;
    Delta_U_int_cond_resid_nJ = 0;
   
end

% Condensation energy at rebound
E_cond_reb_nJ = Delta_U_int_cond_reb_nJ + Delta_E_v_reb_nJ;

% shock wave energy at rebound, Esw,reb
E_SW_reb = U_gas_min - W_visc_max2-E_pot_max2-U_gas_max2-Delta_U_int_cond_reb_nJ/1e9;
E_SW_reb_nJ = E_SW_reb*1e9;
txt1 = ['E_S_W_,_r_e_b = ' num2str((E_SW_reb_nJ)*1e3,'%.2f') ' pJ'];
text(time_Rmax2_us,Work_visc_integral(index_Rmax2)*2e9,txt1,'VerticalAlignment','top');

E_SW_max2_nJ         = E_comp_coll_nJ +E_SW_reb_nJ;   

%% Energy partitioning afterbounce
% Condensation vapor energy for afterbounce
Delta_E_v_resid_nJ = E_v_Rmax2_nJ;
% Condensation energy of afterbounces
E_cond_resid_nJ = Delta_U_int_cond_resid_nJ + Delta_E_v_resid_nJ;

% residual internal energy at Rn3 Ugas,resid
U_gas_resid    = 4*pi/(3*(kappa-1))*((p0+2*sigma/Rn3)*Rn3^3-pg(1)*R0^3);    % residual energy at Rn3
U_gas_resid_nJ = U_gas_resid*1e9;                       % in nJ

% residual potential energy at Rn3 Epot,resid = 0,
% since inner pressure equals external pressure
E_pot_resid_nJ = 0;%4*pi/3*Rn3^3*(p0+2*sigma/Rn3)*1e9;     % in nJ 

% residual energy Eresid
E_resid_nJ = U_gas_resid_nJ + E_pot_resid_nJ;

% viscosity afterbounce
W_visc_afterbounce = Work_visc_integral(end)-Work_visc_integral(index_Rmax2);
W_visc_afterbounce_nJ = W_visc_afterbounce*1e9; % in nJ

% acoustic dissipation afterbounce
W_acoust_afterbounce_nJ = E_B_max2_nJ - E_cond_resid_nJ - E_resid_nJ - W_visc_afterbounce_nJ;

%% Dissipation during entire bubble life
% viscosity
W_visc_total    = Work_visc_integral(end);
W_visc_total_nJ = W_visc_total*1e9; % in nJ

% condensation, U_int_cond_tot + E_cond_tot
U_int_cond_tot_nJ = Delta_U_int_cond_Rmax1_nJ+E_cond_2osc_nJ;
E_cond_total_nJ   = U_int_cond_tot_nJ + E_v_tot_nJ; % in nJ

% shock wave
E_SW_sum_total_nJ = E_SW1_nJ + E_SW_max2_nJ+W_acoust_afterbounce_nJ;% sum of SW1, SW2 and Eacoust

% E_SW_total_nJ = max_E_total_num_nJ-W_visc_total_nJ-E_cond_total_nJ-E_resid_nJ; % to check if E_SW_sum_total_nJ = E_SW_total_nJ

%% Temporal tracking of energy partitioning with Reference level of 100 % as (Eabs - Ev)

E_total_sub_Ev            = U_gas + Work_gas_integral + Delta_U_int_cond; % Reference energy Eabs - Ev in Joule temporal evolution
Uint_div_E_total_sub_Ev   = U_gas./max(E_total_sub_Ev)*100;               % Uint/(Eabs-Ev) in percentage
Delta_U_int_div_E_total_sub_Ev = Delta_U_int_cond./max(E_total_sub_Ev)*100; % Delta_Uint,cond/(Eabs-Ev) in percentage
Wgas_div_E_total_sub_Ev   = Work_gas_integral./max(E_total_sub_Ev)*100;   % Wgas/(Eabs-Ev) in percentage
Wvisc_div_E_total_sub_Ev  = Work_visc_integral./max(E_total_sub_Ev)*100;  % Wvisc/(Eabs-Ev) in percentage
EB_pot_div_E_total_sub_Ev = Work_pot_integral./max(E_total_sub_Ev)*100;   % EB,pot/(Eabs-Ev) in percentage

gcf = figure();
set(gcf,'units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
plot(t*1e6,(Uint_div_E_total_sub_Ev+Delta_U_int_div_E_total_sub_Ev+Wgas_div_E_total_sub_Ev),t*1e6,Uint_div_E_total_sub_Ev,t*1e6,Delta_U_int_div_E_total_sub_Ev,t*1e6,Wgas_div_E_total_sub_Ev,t*1e6,Wvisc_div_E_total_sub_Ev,t*1e6,EB_pot_div_E_total_sub_Ev,'LineWidth',1.5);
xlim([-0.05 1.00]*max(t)*1e6);
ylim([-0.05 1.1]*100);
set(gca,'fontsize',fonts);
xlabel('Time (\mus)');
ylabel('Energy partioning (%)');
legend('E_{abs}-E_{v}','U_{int,cond}','\DeltaU_{int}','W_{gas}','W_{visc}','E_{B,pot}');
title('Energy partitioning with E_{abs}-E_v as 100%');

%%
gcf = figure();
set(gcf,'units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
plot(t*1e6,R*1e6,'LineWidth',1.5)
set(gca,'fontsize',fonts)
xlabel('Time (\mus)')
ylabel('Bubble radius (\mum)')
hold on;
plot(t(1:tind_sw1)*1e6,R(1:tind_sw1)*1e6,'r', t(tind_sw2: tind_sw3)*1e6,R(tind_sw2:tind_sw3)*1e6,'r','LineWidth',2);
hold off;
title('red line marks the R(t) region where p_{in} > p_{out} & acoustics could propogate outward, acoust. simu. should be limited in the marked regions');

%% Save the results in an excel file
% Sheet 1 - Parameters and results
Sheet_Nr_1 = 1;
Sheet1_Cell = cell(22,23);
Sheet1_Cell(1,1:23) = {'pulse duration', num2str(tp),'Bubble expansion','pJ','%','Collapse','pJ','%','Rebound','pJ','%','Afterbounces','pJ','%'...
                      ,'Dissipation during entire bubble life time','pJ','%','Dissipation during bubble expansion','pJ','%'...
                      ,'Dissipation during Bubble collapse and rebound','pJ','%'};  % first line
% Input parameters and output results
Sheet1_Cell(2:31,1) = {'R0(um)','Rnbd/R0','Rncoll1/R0','Rncoll2/R0','tmax1(us)','tmax2(us)','Key results:','Rmax1(um)','Rmin1(um)'...
                      ,'Rmax2(um)','Rmin2(um)','Tosc1(ns)','Tosc2(ns)','Texpansion(ns)','Tcollapse(ns)','Tcoll/Texp','pgas_exp_max(MPa)'...
                      ,'pgas_coll_max(MPa)','U_exp_max(m/s)','U_coll_min(m/s)','U_reb_max(m/s)','R_vap,bd (um)','m_vap,bd (10^-18 kg)'...
                      , 'R_vap,Rmax1 (um)','m_vap,Rmax1 (10^-18 kg)','R_vap,coll1 (um)','m_vap,coll1 (10^-18 kg)','R_vap,Rmax2 (um)'...
                      ,'m_vap,Rmax2 (10^-18 kg)','energy density (kJ/cm^-3)'};%first column
Sheet1_Cell(2,2)    = num2cell(R0*1e6);     % R0(um)
Sheet1_Cell(3,2)    = num2cell(Rn/R0);      % Rnbd/R0
Sheet1_Cell(4,2)    = num2cell(Rn2/R0);     % Rncoll1/R0
Sheet1_Cell(5,2)    = num2cell(Rn3/R0);     % Rncoll2/R0
Sheet1_Cell(6,2)    = num2cell(tRmax*1e6);  % tmax1(us)
Sheet1_Cell(7,2)    = num2cell(tRmax2*1e6); % tmax2(us)
Sheet1_Cell(8,2)    = {''};                 % empty
Sheet1_Cell(9,2)    = num2cell(Rmax*1e6);   % Rmax1(um)
Sheet1_Cell(10,2)   = num2cell(Rmin*1e6);   % Rmin1(um)
Sheet1_Cell(11,2)   = num2cell(Rmax2*1e6);  % Rmax2(um)
Sheet1_Cell(12,2)   = num2cell(Rmin2*1e6);  % Rmin2(um)
Sheet1_Cell(13,2)   = num2cell(Tosc1);      % Tosc1(ns)
Sheet1_Cell(14,2)   = num2cell(Tosc2-Tosc1);% Tosc2(ns)
Sheet1_Cell(15,2)   = num2cell...
                      (time_expansion);     % Texp(ns)
Sheet1_Cell(16,2)   = num2cell...
                      (time_collapse);      % Tcoll(ns) 
Sheet1_Cell(17,2)   = num2cell...
                      (Tcoll_Texp_ratio);   % ratio of tcoll to texpansion
Sheet1_Cell(18,2)   = num2cell...
                      (pgas_exp_max/1e6);   % pgas_exp_max(MPa)
Sheet1_Cell(19,2)   = num2cell...
                      (pgas_coll_max/1e6);  % pgas_coll_max(MPa)
Sheet1_Cell(20,2)   = num2cell...
                      (U_exp_max);          % U_exp_max(m/s)
Sheet1_Cell(21,2)   = num2cell...
                      (U_coll_min);         % U_coll_min(m/s)
Sheet1_Cell(22,2)   = num2cell...
                      (U_reb_max);          % U_reb_max(m/s)
Sheet1_Cell(23,2)   = num2cell...
                      (Rvap_bd*1e6);        % R_vap_bd(um)
Sheet1_Cell(24,2)   = num2cell...
                      (4/3*pi*Rvap_bd^3*rho_vap*1e18); % m_vap_bd(10^-18 kg)
Sheet1_Cell(25,2)   = num2cell...
                      (Rvap_Rmax1*1e6);     % R_vap_Rmax1(um)
Sheet1_Cell(26,2)   = num2cell...
                      (4/3*pi*Rvap_Rmax1^3*rho_vap*1e18); % m_vap_Rmax1(10^-18 kg)   
Sheet1_Cell(27,2)   = num2cell...
                      (Rvap_coll*1e6);      % R_vap_coll(um)
Sheet1_Cell(28,2)   = num2cell...
                      (4/3*pi*Rvap_coll^3*rho_vap*1e18); % m_vap_coll(10^-18 kg)   
Sheet1_Cell(29,2)   = num2cell...
                      (Rvap_Rmax2*1e6);     % R_vap_Rmax2(um)
Sheet1_Cell(30,2)   = num2cell...
                      (4/3*pi*Rvap_Rmax2^3*rho_vap*1e18); % m_vap_Rmax2(10^-18 kg)
Sheet1_Cell(31,2)   = num2cell...
                      (max_E_total_num_nJ*1e-18/(4/3*pi*R0^3)); % energy density (kJ/cm^3)
                  
                  
%% Bubble expansion - energy partitioning: columns C&D&E
Sheet1_Cell(2:16,3) = {'Reference energy  E_abs','Condensation energy E_cond,exp','Condensation vapor energy E_v,exp','Condensation internal energy U_int,cond,exp','Shock wave emission E_SW,bd','Viscous damping W_visc'...
                     ,'Vapor energy at Rmax1, Ev,max1','Internal energy Uint,max1','Potential energy, Epot','Epot from hydrostatic pressure, W_stat','Epot from surface tension,W_surf','Bubble energy at Rmax1, EB,max1', ' ', 'U_int,bd','E_v,bd'}; % column C
Sheet1_Cell(2,4)    = num2cell(max_E_total_num_nJ*1e3);                        %    Reference energy  E_abs in pJ
Sheet1_Cell(2,5)    = num2cell(max_E_total_num_nJ/max_E_total_num_nJ*100);     %    Reference energy  E_abs in %

Sheet1_Cell(3,4)    = num2cell(E_cond_exp_nJ*1e3);                             %    Condensation energy E_cond,exp in pJ
Sheet1_Cell(3,5)    = num2cell(E_cond_exp_nJ/max_E_total_num_nJ*100);          %    Condensation energy E_cond,exp in %

Sheet1_Cell(4,4)    = num2cell(Delta_E_v_exp_nJ*1e3);                          %    Condensation vapor energy Delta_Ev,exp in pJ
Sheet1_Cell(4,5)    = num2cell(Delta_E_v_exp_nJ/max_E_total_num_nJ*100);       %    Condensation vapor energy Delta_Ev,exp in %

Sheet1_Cell(5,4)    = num2cell(Delta_U_int_cond_exp_nJ*1e3);                   %    Condensation internal energy Delta_Uint,cond in pJ
Sheet1_Cell(5,5)    = num2cell(Delta_U_int_cond_exp_nJ/max_E_total_num_nJ*100);%    Condensation internal energy Delta_Uint,cond in %

Sheet1_Cell(6,4)    = num2cell(E_SW1_nJ*1e3);                                  %    Shock wave emission E_SW,max1 in pJ
Sheet1_Cell(6,5)    = num2cell(E_SW1_nJ/max_E_total_num_nJ*100);               %    Shock wave emission E_SW,max1 in %

Sheet1_Cell(7,4)    = num2cell(max_Work_visc_integral_nJ*1e3);                 %    Viscous damping W_visc in pJ
Sheet1_Cell(7,5)    = num2cell(max_Work_visc_integral_nJ...
                               /max_E_total_num_nJ*100);                       %    Viscous damping W_visc in %

Sheet1_Cell(8,4)    = num2cell(E_v_Rmax1_nJ*1e3);                              %    Vapor energy at Rmax1, Ev,Rmax1 in pJ
Sheet1_Cell(8,5)    = num2cell(E_v_Rmax1_nJ/max_E_total_num_nJ*100);           %    Vapor energy at Rmax1, Ev,Rmax1 in %

Sheet1_Cell(9,4)    = num2cell(U_gas_max1_nJ*1e3);                             %    Internal energy Uint,max1 in pJ
Sheet1_Cell(9,5)    = num2cell(U_gas_max1_nJ/max_E_total_num_nJ*100);          %    Internal energy Uint,max1 in %

Sheet1_Cell(10,4)   = num2cell(max_Work_pot_integral_nJ*1e3);                  %    Epot from hydrostatic pressure, W_stat in pJ
Sheet1_Cell(10,5)   = num2cell(max_Work_pot_integral_nJ...
                               /max_E_total_num_nJ*100);                       %    Epot from hydrostatic pressure, W_stat in %

Sheet1_Cell(11,4)   = num2cell(max_Work_stat_nJ*1e3);                          %    Epot from hydrostatic pressure, W_stat in pJ
Sheet1_Cell(11,5)   = num2cell(max_Work_stat_nJ/max_E_total_num_nJ*100);       %    Epot from hydrostatic pressure, W_stat in %

Sheet1_Cell(12,4)   = num2cell(max_Work_surf_integral_nJ*1e3);                 %    Epot from surface tension,W_surf in pJ
Sheet1_Cell(12,5)   = num2cell(max_Work_surf_integral_nJ...
                               /max_E_total_num_nJ*100);                       %    Epot from surface tension,W_surf in %


Sheet1_Cell(13,4)   = num2cell(E_B_max1_nJ*1e3);                               %    Bubble energy at Rmax1, EB,max1 in pJ
Sheet1_Cell(13,5)   = num2cell(E_B_max1_nJ/max_E_total_num_nJ*100);            %    Bubble energy at Rmax1, EB,max1 in %

Sheet1_Cell(15,4)   = num2cell(U_int_bd_pJ);                                   %    Internal energy after breakdown in pJ
Sheet1_Cell(15,5)   = num2cell(U_int_bd_per);                                  %    Internal energy after breakdown in %

Sheet1_Cell(16,4)   = num2cell(E_v_bd_pJ);                                     %    Vaporization energy after breakdown in pJ
Sheet1_Cell(16,5)   = num2cell(E_v_bd_per);                                    %    Vaporization energy after breakdown in %

%% Collapse         - energy partitioning: columns F,G,H
Sheet1_Cell(2:10,6) = {'Reference energy  E_B_max1','Condensation energy E_cond,coll','Condensation vapor energy E_v,coll','Condensation internal energy U_int,cond,coll','Viscous damping W_visc','Vapor energy at Rmin, Ev,min1'...
                     ,'Internal energy of bubble, U_int,min','Internal energy of liquid, E_compr,coll','Compression energy at Rmin1, E_compr,total'}; % column F

Sheet1_Cell(2,7)    = num2cell(E_B_max1_nJ*1e3);                               %    Reference energy  E_B_max1 in pJ
Sheet1_Cell(2,8)    = num2cell(E_B_max1_nJ...
                               /E_B_max1_nJ*100);                              %    Reference energy  E_B_max1 in %                 

Sheet1_Cell(3,7)    = num2cell(E_cond_coll_nJ*1e3);                            %    Condensation energy E_cond,coll in pJ
Sheet1_Cell(3,8)    = num2cell(E_cond_coll_nJ/E_B_max1_nJ*100);                %    Condensation energy E_cond,coll in %

Sheet1_Cell(4,7)    = num2cell(Delta_E_v_coll_nJ*1e3);                         %    Condensation vapor energy Delta_E_v_coll in pJ
Sheet1_Cell(4,8)    = num2cell(Delta_E_v_coll_nJ/E_B_max1_nJ*100);             %    Condensation vapor energy Delta_E_v_coll in %

Sheet1_Cell(5,7)    = num2cell(Delta_U_int_cond_coll_nJ*1e3);                  %    Condensation internal energy U_int,cond,coll in pJ
Sheet1_Cell(5,8)    = num2cell(Delta_U_int_cond_coll_nJ/E_B_max1_nJ*100);      %    Condensation internal energy U_int,cond,coll in %

Sheet1_Cell(6,7)    = num2cell(W_visc_min_nJ*1e3);                             %    Viscous damping W_visc in pJ
Sheet1_Cell(6,8)    = num2cell(W_visc_min_nJ/E_B_max1_nJ*100);                 %    Viscous damping W_visc in %
                            
Sheet1_Cell(7,7)    = num2cell(E_v_coll_nJ*1e3);                               %    Vapor energy at Rmin Ev,min1 in pJ
Sheet1_Cell(7,8)    = num2cell(E_v_coll_nJ/E_B_max1_nJ*100);                   %    Vapor energy at Rmin Ev,min1 in %

Sheet1_Cell(8,7)    = num2cell(U_gas_min_nJ*1e3);                              %    Internal energy of bubble, U_int,min in pJ
Sheet1_Cell(8,8)    = num2cell(U_gas_min_nJ/E_B_max1_nJ*100);                  %    Internal energy of bubble, U_int,min in %

Sheet1_Cell(9,7)    = num2cell(E_comp_coll_nJ*1e3);                            %    Internal energy of liquid, E_compr,coll in pJ
Sheet1_Cell(9,8)    = num2cell(E_comp_coll_nJ/E_B_max1_nJ*100);                %    Internal energy of liquid, E_compr,coll in %

E_compr_total_nJ    = E_comp_coll_nJ+U_gas_min_nJ+E_v_coll_nJ;                 %    Compression energy at Rmin1, E_compr,total in nJ
Sheet1_Cell(10,7)   = num2cell((E_compr_total_nJ)*1e3);                        %    Compression energy at Rmin1, E_compr,total in pJ
Sheet1_Cell(10,8)   = num2cell((E_compr_total_nJ)...
                                /E_B_max1_nJ*100);                             %    Compression energy at Rmin1, E_compr,total in %

%% Rebound          - energy partitioning: columns I J K
Sheet1_Cell(2:15,9) = {'Reference energy E_compr,total','Condensation energy E_cond,reb','Condensation vapor energy E_v,reb','Condensation internal energy U_int,cond,reb','Shock wave emission E_SW,reb','Shock wave from compressed liquid E_SWL',...
                       'Shock wave from bubble interior E_SWB','Viscous damping W_visc','Vapor energy at Rmax2, Ev,max2',...
                       'Internal energy at Rmax2, U_int,max2','Potential energy, Epot2','Epot from hydrostatic pressure, W_stat','Epot from surface tension, W_surf',...
                       'Bubble energy at Rmax2, E_B,max2'};                     %    column I
Sheet1_Cell(2,10)    = num2cell(E_compr_total_nJ*1e3);                          %    Reference energy E_compr,total in pJ
Sheet1_Cell(2,11)    = num2cell(E_compr_total_nJ/...
                                E_compr_total_nJ*100);                          %    Reference energy E_compr,total in %

Sheet1_Cell(3,10)     = num2cell(E_cond_reb_nJ*1e3);                             %    Condensation energy E_cond,reb in pJ
Sheet1_Cell(3,11)     = num2cell(E_cond_reb_nJ/E_compr_total_nJ*100);          %    Condensation energy E_cond,reb in %

Sheet1_Cell(4,10)     = num2cell(Delta_E_v_reb_nJ*1e3);                          %    Condensation vapor energy Delta_Ev,reb in pJ
Sheet1_Cell(4,11)     = num2cell(Delta_E_v_reb_nJ/E_compr_total_nJ*100);       %    Condensation vapor energy Delta_Ev,reb in %
                            
                            
Sheet1_Cell(5,10)    = num2cell(Delta_U_int_cond_reb_nJ*1e3);                   %    Condensation internal energy U_int,cond,reb in pJ
Sheet1_Cell(5,11)    = num2cell(Delta_U_int_cond_reb_nJ/E_compr_total_nJ*100);  %    Condensation internal energy U_int,cond,reb in %                            
                            
Sheet1_Cell(6,10)    = num2cell(E_SW_max2_nJ*1e3);                              %    Shock wave emission E_SW,max2 in pJ
Sheet1_Cell(6,11)    = num2cell(E_SW_max2_nJ/E_compr_total_nJ*100);             %    Shock wave emission E_SW,max2 in %
Sheet1_Cell(7,10)    = num2cell(E_comp_coll_nJ*1e3);                            %    Shock wave from compressed liquid E_compr,coll in pJ
Sheet1_Cell(7,11)    = num2cell(E_comp_coll_nJ/E_compr_total_nJ*100);           %    Shock wave from compressed liquid E_compr,coll in %
Sheet1_Cell(8,10)    = num2cell(E_SW_reb_nJ*1e3);                               %    Shock wave from bubble interior E_SW,B in pJ
Sheet1_Cell(8,11)    = num2cell(E_SW_reb_nJ/E_compr_total_nJ*100);              %    Shock wave from bubble interior E_SW,B in %

Sheet1_Cell(9,10)    = num2cell(W_visc_max2_nJ*1e3);                            %    Viscous damping W_visc in pJ
Sheet1_Cell(9,11)    = num2cell(W_visc_max2_nJ/E_compr_total_nJ*100);           %    Viscous damping W_visc in %

Sheet1_Cell(10,10)   = num2cell(E_v_Rmax2_nJ*1e3);                              %    Vapor energy at Rmax2 E_v_Rmax2 in pJ
Sheet1_Cell(10,11)   = num2cell(E_v_Rmax2_nJ/E_compr_total_nJ*100);             %    Vapor energy at Rmax2 E_v_Rmax2 in %

Sheet1_Cell(11,10)   = num2cell(U_gas_max2_nJ*1e3);                             %    Internal energy U_int,max2 in pJ
Sheet1_Cell(11,11)   = num2cell(U_gas_max2_nJ/E_compr_total_nJ*100);            %    Internal energy U_int,max2 in %

Sheet1_Cell(12,10)    = num2cell(E_pot_max2_nJ*1e3);                             %    Epot2 in pJ
Sheet1_Cell(12,11)    = num2cell(E_pot_max2_nJ/E_compr_total_nJ*100);            %    Epot2 in %
Sheet1_Cell(13,10)    = num2cell(W_stat_max2_nJ*1e3);                            %    Epot from hydrostatic pressure, W_stat in pJ
Sheet1_Cell(13,11)    = num2cell(W_stat_max2_nJ/E_compr_total_nJ*100);           %    Epot from hydrostatic pressure, W_stat in %
Sheet1_Cell(14,10)    = num2cell(W_surf_max2_nJ*1e3);                            %    Epot from surface tension, W_surf in pJ
Sheet1_Cell(14,11)    = num2cell(W_surf_max2_nJ/E_compr_total_nJ*100);           %    Epot from surface tension, W_surf in %

Sheet1_Cell(15,10)    = num2cell(E_B_max2_nJ*1e3);                               %    Bubble energy at Rmax2, E_B,max2 in pJ
Sheet1_Cell(15,11)    = num2cell(E_B_max2_nJ/E_compr_total_nJ*100);              %    Bubble energy at Rmax2, E_B,max2 in %

%% Afterbound/residual bubble - energy partitioning: columns L,M,N
Sheet1_Cell(2:10,12) = {'Reference energy E_B,max2','Condensation energy E_cond,resid','Condensation vapor energy E_v,resid','Condensation internal energy U_int,cond,resid',...
                       'Acoustic emission E_acoust','Viscous damping W_visc','Residual bubble energy E_resid','Residual bubble internal energy U_resid',...
                       'Residual bubble potential energy E_pot_resid'}; % column L
Sheet1_Cell(2,13)    = num2cell(E_B_max2_nJ*1e3);                               %    Reference energy E_B,max2 in pJ
Sheet1_Cell(2,14)    = num2cell(E_B_max2_nJ/...
                                E_B_max2_nJ*100);                               %    Reference energy E_B,max2 in %

Sheet1_Cell(3,13)    = num2cell(E_cond_resid_nJ*1e3);                           %    Condensation energy E_cond,resid in pJ
Sheet1_Cell(3,14)    = num2cell(E_cond_resid_nJ/...
                                E_B_max2_nJ*100);                               %    Condensation energy E_cond,resid in %                            

Sheet1_Cell(4,13)    = num2cell(Delta_E_v_resid_nJ*1e3);                        %    Condensation vapor energy E_v,resid in pJ
Sheet1_Cell(4,14)    = num2cell(Delta_E_v_resid_nJ/...
                                E_B_max2_nJ*100);                               %    Condensation vapor energy E_v,resid in %                            

Sheet1_Cell(5,13)    = num2cell(Delta_U_int_cond_resid_nJ*1e3);                 %    Condensation internal energy U_int,cond,resid in pJ
Sheet1_Cell(5,14)    = num2cell(Delta_U_int_cond_resid_nJ/...
                                E_B_max2_nJ*100);                               %    Condensation internal energy U_int,cond,resid in % 
                          
Sheet1_Cell(6,13)    = num2cell(W_acoust_afterbounce_nJ*1e3);                   %    Acoustic emission E_acoust in pJ
Sheet1_Cell(6,14)    = num2cell(W_acoust_afterbounce_nJ/...
                                E_B_max2_nJ*100);                               %    Acoustic emission E_acoust in %
                            
Sheet1_Cell(7,13)    = num2cell(W_visc_afterbounce_nJ*1e3);                     %    Viscous damping W_visc in pJ
Sheet1_Cell(7,14)    = num2cell(W_visc_afterbounce_nJ/...
                                E_B_max2_nJ*100);                               %    Viscous damping W_visc in %
                            
Sheet1_Cell(8,13)    = num2cell(E_resid_nJ*1e3);                                %    Residual bubble energy E_resid in pJ
Sheet1_Cell(8,14)    = num2cell(E_resid_nJ/...
                                E_B_max2_nJ*100);                               %    Residual bubble energy E_resid in % 
                            
Sheet1_Cell(9,13)    = num2cell(U_gas_resid_nJ*1e3);                            %    Residual bubble internal energy U_resid in pJ
Sheet1_Cell(9,14)    = num2cell(U_gas_resid_nJ/...
                                E_B_max2_nJ*100);                               %    Residual bubble internal energy U_resid in % 

Sheet1_Cell(10,13)    = num2cell(E_pot_resid_nJ*1e3);                           %    Residual bubble potential energy E_pot_resid in pJ
Sheet1_Cell(10,14)    = num2cell(E_pot_resid_nJ/...
                                E_B_max2_nJ*100);                               %    Residual bubble potential energy E_pot_resid in %                             

%% Dissipation during entire bubble life time: column O,P,Q
Sheet1_Cell(2:12,15) = {'Reference energy Eabs','Shock wave emission ESW+Eacoust','Condensation Econd',...
                       'Viscous damping Wvisc','Residual bubble energy E_residual',' ', ...
                       'Vaporization and condensation', 'E_v,bd', 'Delta_E_v_exp', 'Delta_E_v_coll&reb','E_v,Rmax2'}; % column O
Sheet1_Cell(2,16)    = num2cell(max_E_total_num_nJ*1e3);                        %    Reference energy Eabs in pJ
Sheet1_Cell(2,17)    = num2cell(max_E_total_num_nJ/...
                                max_E_total_num_nJ*100);                        %    Reference energy Eabs in %    
                            
Sheet1_Cell(3,16)    = num2cell(E_SW_sum_total_nJ*1e3);                         %    Shock wave emission ESW+Eacoust in pJ
Sheet1_Cell(3,17)    = num2cell(E_SW_sum_total_nJ/...
                                max_E_total_num_nJ*100);                        %    Shock wave emission ESW+Eacoust in %  
                            
Sheet1_Cell(4,16)    = num2cell(E_cond_total_nJ*1e3);                           %    Condensation total Econd in pJ
Sheet1_Cell(4,17)    = num2cell(E_cond_total_nJ/...
                                max_E_total_num_nJ*100);                        %    Condensation total Econd in %
                            
Sheet1_Cell(5,16)    = num2cell(W_visc_total_nJ*1e3);                           %    Viscous damping Wvisc in pJ
Sheet1_Cell(5,17)    = num2cell(W_visc_total_nJ/...
                                max_E_total_num_nJ*100);                        %    Viscous damping Wvisc in % 
                            
Sheet1_Cell(6,16)    = num2cell(E_resid_nJ*1e3);                                %    Residual bubble energy E_resid in pJ
Sheet1_Cell(6,17)    = num2cell(E_resid_nJ/...
                                max_E_total_num_nJ*100);                        %    Residual bubble energy E_resid in % 
 
Sheet1_Cell(9,16)    = num2cell(E_v_tot*1e12);                                  %    E_v,bd in pJ
Sheet1_Cell(9,17)    = num2cell(E_v_tot/E_v_tot*100);                           %    E_v,bd in %  

Sheet1_Cell(10,16)    = num2cell(Delta_E_v_exp*1e12);                           %    Delta_E_v,exp in pJ
Sheet1_Cell(10,17)    = num2cell(Delta_E_v_exp/E_v_tot*100);                    %    Delta_E_v,exp in % 

Sheet1_Cell(11,16)    = num2cell(Delta_E_v_coll_reb*1e12);                      %    Delta_E_v,coll&reb in pJ
Sheet1_Cell(11,17)    = num2cell(Delta_E_v_coll_reb/E_v_tot*100);               %    Delta_E_v,coll&reb in % 

Sheet1_Cell(12,16)    = num2cell(E_v_Rmax2*1e12);                               %    Ev,Rmax2 in pJ
Sheet1_Cell(12,17)    = num2cell(E_v_Rmax2/E_v_tot*100);                        %    Ev,Rmax2 in %
                            
%% Dissipation during bubble expansion: column R,S,T
Sheet1_Cell(2:6,18) = {'Reference energy Eabs','Shock wave emission ESW,1','Condensation Econd,exp',...
                       'Viscous damping Wvisc','Bubble energy at Rmax1, EB,max1'}; % column R
Sheet1_Cell(2,19)    = num2cell(max_E_total_num_nJ*1e3);                        %    Reference energy Eabs in pJ
Sheet1_Cell(2,20)    = num2cell(max_E_total_num_nJ/...
                                max_E_total_num_nJ*100);                        %    Reference energy Eabs in %
                            
Sheet1_Cell(3,19)    = num2cell(E_SW1_nJ*1e3);                                  %    Shock wave emission ESW,1 in pJ
Sheet1_Cell(3,20)    = num2cell(E_SW1_nJ/...
                                max_E_total_num_nJ*100);                        %    Shock wave emission ESW,1 in %
                                                  
Sheet1_Cell(4,19)    = num2cell(E_cond_exp_nJ*1e3);                             %    Condensation during expansion Econd,exp in pJ
Sheet1_Cell(4,20)    = num2cell(E_cond_exp_nJ/...
                                max_E_total_num_nJ*100);                        %    Condensation during expansion Econd,exp in %
 
Sheet1_Cell(5,19)    = num2cell(max_Work_visc_integral_nJ*1e3);                 %    Viscous damping Wvisc in pJ
Sheet1_Cell(5,20)    = num2cell(max_Work_visc_integral_nJ/...
                                max_E_total_num_nJ*100);                        %    Viscous damping Wvisc in %
                                
Sheet1_Cell(6,19)    = num2cell(E_B_max1_nJ*1e3);                               %    Bubble energy at Rmax1, EB,max1 in pJ
Sheet1_Cell(6,20)    = num2cell(E_B_max1_nJ/...
                                max_E_total_num_nJ*100);                        %    Bubble energy at Rmax1, EB,max1 in % 

%% Dissipation during Bubble collapse and rebound: columns U,V,W
Sheet1_Cell(2:6,21) = {'Reference energy EB,max1','Shock wave emission ESW,2','Condensation Econd,coll + Econd,reb',...
                       'Viscous damping Wvisc','Bubble energy at Rmax2, EB,max2'}; % column R
Sheet1_Cell(2,22)    = num2cell(E_B_max1_nJ*1e3);                               %    Reference energy EB,max1 in pJ
Sheet1_Cell(2,23)    = num2cell(E_B_max1_nJ/...
                                E_B_max1_nJ*100);                               %    Reference energy EB,max1 in %                    

Sheet1_Cell(3,22)    = num2cell(E_SW_max2_nJ*1e3);                              %    Shock wave emission ESW,2 in pJ
Sheet1_Cell(3,23)    = num2cell(E_SW_max2_nJ/...
                                E_B_max1_nJ*100);                               %    Shock wave emission ESW,2 in %   

Sheet1_Cell(4,22)    = num2cell((E_cond_coll_nJ+E_cond_reb_nJ)*1e3);            %    Condensation & heat conduction Econd in pJ
Sheet1_Cell(4,23)    = num2cell((E_cond_coll_nJ+E_cond_reb_nJ)/...
                                E_B_max1_nJ*100);                               %    Condensation & heat conduction Econd in %
                            
Sheet1_Cell(5,22)    = num2cell((W_visc_min_nJ+W_visc_max2_nJ)*1e3);            %    Viscous damping Wvisc in pJ
Sheet1_Cell(5,23)    = num2cell((W_visc_min_nJ+W_visc_max2_nJ)/...
                                E_B_max1_nJ*100);                               %    Viscous damping Wvisc in % 
                            
Sheet1_Cell(6,22)    = num2cell(E_B_max2_nJ*1e3);                               %    Bubble energy at Rmax2, EB,max2 in pJ
Sheet1_Cell(6,23)    = num2cell(E_B_max2_nJ/...
                                E_B_max1_nJ*100);                        %    Bubble energy at Rmax2, EB,max2 in % 

% Compatibility issue to write a cell in excile file
if verLessThan('matlab','9.6')                      
    xlswrite(filename,Sheet1_Cell,Sheet_Nr_1,'A1');      % xlswrite for matlab before R2019a (Matlab 9.6)
else
    writecell(Sheet1_Cell, filename,'Sheet',Sheet_Nr_1); % writecell is introduced in R2019a (Matlab 9.6)
end

% Sheet 2 - Raw data
Sheet2_Cell = cell(length(t)+1,13);
Sheet2_Cell(1,:) = {'time(us)', 'radius(um)','Rnt(um)','Rvdw(um)','U(m/s)','p_gas(MPa)','E_total(nJ)','U_gas(nJ)','W_gas(nJ)','W_visc(nJ)','W_surf(nJ)','E_pot(nJ)','Empty'};
Sheet2_Cell(2:end,1) = num2cell(t*1e6);   % time in us
Sheet2_Cell(2:end,2) = num2cell(R*1e6);   % radius in um
Sheet2_Cell(2:end,3) = num2cell(Rnt*1e6); % Rn in um
Sheet2_Cell(2:end,4) = num2cell(Rvdw*1e6);% Rvdw in um
Sheet2_Cell(2:end,5) = num2cell(U);       % bubble wall velocity U in m/s
Sheet2_Cell(2:end,6) = num2cell(pg/1e6);  % gas pressure in MPa
Sheet2_Cell(2:end,7) = num2cell(E_total_num*1e9); % total energy in nJ
Sheet2_Cell(2:end,8) = num2cell(U_gas*1e9);       % internal energy in nJ
Sheet2_Cell(2:end,9) = num2cell(Work_gas_integral*1e9);   % work done by the gas in nJ
Sheet2_Cell(2:end,10) = num2cell(Work_visc_integral*1e9); % work done by viscocity in nJ
Sheet2_Cell(2:end,11) = num2cell(Work_surf_integral*1e9); % work done by surface tension in nJ
Sheet2_Cell(2:end,12) = num2cell(Work_pot_integral*1e9);  % potential energy in nJ
Sheet2_Cell(2:end,13) = num2cell(0);      % empty
Sheet_Nr_2 = 2;

% Compatibility issue to write a cell in excile file, the same as with Sheet 1
if verLessThan('matlab','9.6')
    xlswrite(filename,Sheet2_Cell,Sheet_Nr_2,'A1');       % before R2019a
else
    writecell(Sheet2_Cell, filename,'Sheet',Sheet_Nr_2);  % R2019a and after
end


else
    disp('---            Energy partitioning is switched off           ---');
end

%% Assign output
varargout{1} = t;
varargout{2} = R;
varargout{3} = pg;
varargout{4} = U;
varargout{5} = tind_sw;
varargout{6} = index_Rmin;
end

