function [varargout] = SelcurvReshaping(t,tindstart,tindend,dtkb,r_container,u_container,p_container,Rtt,Utt,Ptt,varargin)
% This function simulates acoustic radiation in the liquid around the
% bubble using the Kirkwood-Bethe hyposis and reshapes the multivalued shock wave front, See Ref:
% Liang et al. J. Fluid. Mech. 940, A5 (2022). DOI: 10.1017/jfm.2022.202

% Usage:
% [varargout] =
% SelcurvReshaping(t,tindstart,tindend,dtkb,r_container,u_container,p_container,Rtt,Utt,Ptt)
%
% Explicit input arguments [in SI units]:
% t     - time axis for bubble dynamics
% tindstart - index of time array t() where simulation of sw starts
% tindend   - index of time array t() where simulation of sw ends
% dtkb      - time step for the acoustic radiation using the KB algorithm
% r_container - propagation distance, array size [T_KB,T_Gilmore]
% u_container - local velocity along propagating distance, array size [T_KB,T_Gilmore]
% p_container - pressure along propagating distance, array size [T_KB,T_Gilmore]
% Rtt - bubble radius interpolated on the KB time axis
% Utt - bubble wall velocity interpolated on the KB time axis
% Ptt - pressure inside the bubble interpolated on the KB time axis
%
% Optional parameters
% time_q - queried time instants in ns
%
% Output arguments:
% varargout{1} - original curves
% varargout{2} - reshaped curves, reshaping p(r) from reshaped u(r)
% varargout{3} - reshaped curves, reshaping p(r) and u(r) seperately
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

defaultnquery     = 1;
defaulttime_q     = t(round((tindstart+tindend)/2))*1e9; % in ns
defaultisExtrafig = 0;  % plot extra figure

validScalarPosNum    = @(x) isnumeric(x) && isscalar(x) && (x > 0);

addParameter (p,'nquery',defaultnquery,validScalarPosNum);
addParameter (p,'time_q',defaulttime_q);
addParameter (p,'isExtrafig', defaultisExtrafig);

parse(p,varargin{:});   % Parse param.-value pairs
param = p.Results;      % transfer res. to structure
clear p                 % delete the object p that is not needed


time_q = param.time_q;  % querried time instants in unit ns
nquery = length(time_q);% number of queried curves

isExtrafig = param.isExtrafig; % extra figure for reshaping p(r) and u(r) curves respectively

assert(min(time_q)>= t(tindstart)*1e9 && max(time_q)<=t(tindend)*1e9,'error in querried time, check your input time_q'); % querried time must be limited in the time window t(tindstart:tindend)

disp('---      reshaping shock front at selected time instants     ---');

fontsize = 12; 
%% time index transformation
% time index of the querried time in array T_Gilmore
tindg = interp1(t(tindstart:tindend)*1e9, tindstart:tindend, time_q,'nearest');

% time index of the querried time in the acoustic emission time axis T_kb
ttind = round((t(round(tindg))-t(tindstart))/dtkb)+1; 

%% preallocate arrays
rsize = size(r_container);                              % the size of the array (T_KB,T_Gilmore)
p_u_r_n_curves            = NaN(rsize(1,2),(nquery*3)); % original curves
p_r_u_r_n_curves_reshaped = NaN(rsize(1,2),(nquery*4)); % reshaping p_r and u_r respectively
p_r_von_u_r_reshaped      = NaN(size(p_u_r_n_curves));  % reshaping p_r according to reshaped u_r curves

%% original u(r) and p(r) curves
figure()
hold on
for n = 1:nquery
    rplot = r_container(ttind(n),:)*1e6; % [1 x T_Gilmore] um
    pplot = p_container(ttind(n),:)/1e6; % [1 x T_Gilmore] MPa
    uplot = u_container(ttind(n),:);     % [1 x T_Gilmore] m/s
    rplot(isnan(rplot)) = 0;    % we must get rid of the NaN's
    pplot(isnan(pplot)) = 0;    % we must get rid of the NaN's
    uplot(isnan(uplot)) = 0;    % we must get rid of the NaN's
    j = find(rplot);            % and only find the valid values (nonzeros)
    
    txt_time_info = ['p(r)@t = ' num2str(time_q(n),'%.2f'),'ns']; % for displaying in legend
    yyaxis left
    plot([rplot(j) Rtt(ttind(n))*1e6], [pplot(j) Ptt(ttind(n))/1e6],'LineWidth',1.5,'DisplayName',txt_time_info)     % Rtt,Ptt and Utt are the values at bubble wall
    xlabel('r(\mum)','fontsize',fontsize)
    ylabel('p(MPa)','fontsize',fontsize)   

    txt_time_info = ['u(r)@t = ' num2str(time_q(n),'%.2f'),'ns'];
    yyaxis right
    plot([rplot(j) Rtt(ttind(n))*1e6], [uplot(j) Utt(ttind(n))],'LineWidth',1.5,'DisplayName',txt_time_info)
    ylabel('u (m/s))','fontsize',fontsize)
    set(gca,'fontsize',fontsize)
    title('original curves')

    % fliping the array up side down to display progagating distance from
    % small to large, see the contour map of r(t_Gilmore,t_KB) in Figure 12 
    rplot_temp_um  = flipud([rplot(j) Rtt(ttind(n))*1e6]'); % [T_Gilmore x 1] 
    pplot_temp_MPa = flipud([pplot(j) Ptt(ttind(n))/1e6]'); % [T_Gilmore x 1]
    uplot_temp     = flipud([uplot(j) Utt(ttind(n))]');     % [T_Gilmore x 1]

    % save results for the original curves
    p_u_r_n_curves(1:length(rplot_temp_um), 3*n-2)  = rplot_temp_um;
    p_u_r_n_curves(1:length(pplot_temp_MPa),3*n-1)  = pplot_temp_MPa;
    p_u_r_n_curves(1:length(uplot_temp),    3*n)    = uplot_temp;
    
end
hold off
legend show

%% reshaping p(r) and u(r) curves respectively
if isExtrafig == 1
figure ()
hold on
for n = 1:nquery
    r_temp = p_u_r_n_curves(:,3*n-2); % original curves
    p_temp = p_u_r_n_curves(:,3*n-1);
    u_temp = p_u_r_n_curves(:,3*n);
    
    [r_p_reshaped,p_reshaped,~,~] = shockFrontreshaping(r_temp,p_temp); % reshaping p(r) curves
    [r_u_reshaped,u_reshaped,~,~] = shockFrontreshaping(r_temp,u_temp); % reshaping u(r) curves
    
    p_r_u_r_n_curves_reshaped(1:length(r_p_reshaped),4*n-3)    = r_p_reshaped; % save results
    p_r_u_r_n_curves_reshaped(1:length(p_reshaped),  4*n-2)    = p_reshaped;
    p_r_u_r_n_curves_reshaped(1:length(r_u_reshaped),4*n-1)    = r_u_reshaped;
    p_r_u_r_n_curves_reshaped(1:length(u_reshaped),  4*n)      = u_reshaped;
    
    txt_time_info = ['p(r)@t = ' num2str(time_q(n),'%.2f'),'ns'];
    
    yyaxis left
    plot(r_p_reshaped, p_reshaped,'LineWidth',1.5,'DisplayName',txt_time_info)
    xlabel('r(\mum)','fontsize',fontsize)
    ylabel('p(MPa)','fontsize',fontsize)
    
    txt_time_info = ['u(r)@t = ' num2str(time_q(n),'%.2f'),'ns'];
    
    yyaxis right
    plot(r_u_reshaped, u_reshaped,'LineWidth',1.5,'DisplayName',txt_time_info)
    ylabel('u (m/s))','fontsize',fontsize)
    set(gca,'fontsize',fontsize)
    title('reshaping u(r) and p(r) respectively')
         
end
hold off
legend show
end
%% reshaping p(r) curves according to the reshaped u(r) curves
figure ()
hold on
for n = 1:nquery
    
r_temp = p_u_r_n_curves(:,3*n-2); % original curves
p_temp = p_u_r_n_curves(:,3*n-1);
u_temp = p_u_r_n_curves(:,3*n);

% reshaping p(r) curves according to the reshaped u(r) curves
[r_reshaped,u_reshaped,p_reshaped,~,~] = shockFronttransfer(r_temp,u_temp,p_temp); % rout ~ and pout ~ for the peak of the shock front

p_r_von_u_r_reshaped(1:length(r_reshaped),3*n-2)    = r_reshaped; % save results
p_r_von_u_r_reshaped(1:length(u_reshaped),3*n-1)    = u_reshaped;
p_r_von_u_r_reshaped(1:length(p_reshaped),3*n)      = p_reshaped;

    txt_time_info = ['p(r)@t = ' num2str(time_q(n),'%.2f'),'ns'];
    yyaxis left
    plot(r_reshaped, p_reshaped,'LineWidth',1.5,'DisplayName',txt_time_info)
    xlabel('r(\mum)','fontsize',fontsize)
    ylabel('p(MPa)','fontsize',fontsize)
    
    txt_time_info = ['u(r)@t = ' num2str(time_q(n),'%.2f'),'ns'];
    yyaxis right
    plot(r_reshaped, u_reshaped,'LineWidth',1.5,'DisplayName',txt_time_info)
    ylabel('u (m/s))','fontsize',fontsize)
    set(gca,'fontsize',fontsize)
    title('reshaped curves')
    
end
hold off
legend show

varargout{1} = p_u_r_n_curves;            % original p(r) and u(r) curves
varargout{2} = p_r_von_u_r_reshaped;      % reshaping p(r) curves according to the reshaped u(r) curves
varargout{3} = p_r_u_r_n_curves_reshaped; % reshaping p(r) and u(r) curves seperately
end