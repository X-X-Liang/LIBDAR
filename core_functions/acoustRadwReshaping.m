function [varargout] = acoustRadwReshaping(t,R,tp,tindstart,tindend,tindRmin,dtkb,r_container,u_container,p_container,Rtt,Utt,Ptt,varargin)
% This function simulates acoustic radiation in the liquid around the bubble using the Kirkwood-Bethe hyposis 
% reshapes the multivalued shock wave front, and creates animation and plotting. See Ref:
% Liang et al. J. Fluid. Mech. 940, A5 (2022). DOI: 10.1017/jfm.2022.202

% Usage:
% [varargout] =
% acoustRadwReshaping(t,tindstart,tindend,tindRmin,dtkb,r_container,u_container,p_container,Rtt,Utt,Ptt)
%
% Explicit input arguments [in SI units]:
% t     - time axis for bubble dynamics
% R     - bubble radius
% tp    - pulse duration
% tindstart - index of time array t() where simulation of sw starts
% tindend   - index of time array t() where simulation of sw ends
% tindRmin  - index of time array t() where Rmin is reached
% dtkb      - time step for the acoustic radiation using the KB algorithm
% r_container - propagation distance, array size [T_KB,T_Gilmore]
% u_container - local velocity along propagating distance, array size [T_KB,T_Gilmore]
% p_container - pressure along propagating distance, array size [T_KB,T_Gilmore]
% Rtt - bubble radius interpolated on the KB time axis
% Utt - bubble wall velocity interpolated on the KB time axis
% Ptt - pressure inside the bubble interpolated on the KB time axis
%
% Optional parameters
% Nsel    - nr. of curve interval for visulization
% isARexp - is acousti radiation at expansion phase
% xylim   - x and y limits for visulization
% isExtrafig - is extra figure needed
% time_q     - queried time instants
%
% Output arguments:
% varargout{1} - 
% varargout{2} - 
% varargout{3} - 
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

defaultNsel       = 4; % select every 4 curves for plotting and animation
defaultisARexp    = 1; % expansion phase
defaultxylim      = [1 100 0.1 1e4]; % xmin = 1um, xmax=100 um, ymin = 0.1MPa, ymax = 1e4 MPa
defaultisExtrafig = 0; % do not plot extra figure
defaulttime_q     = t(round((tindstart+tindend)/2))*1e9; % in unit ns

validScalarPosNum    = @(x) isnumeric(x) && isscalar(x) && (x > 0);

addParameter (p,'Nsel',defaultNsel,validScalarPosNum);
addParameter (p,'isARexp',defaultisARexp);
addParameter (p,'xylim',defaultxylim);
addParameter (p,'isExtrafig', defaultisExtrafig);
addParameter (p,'time_q',defaulttime_q)

parse(p,varargin{:});   % Parse param.-value pairs
param = p.Results;      % transfer res. to structure
clear p                 % delete the object p that is not needed

Nsel    = param.Nsel;
isARexp = param.isARexp;
xylim   = param.xylim;
isExtrafig = param.isExtrafig;
time_q     = param.time_q;  % querried time instants in unit ns


%% constants
fontsize = 22; 

%% preallocate to store results for animation Res_R_um_t_us, Res_P_MPa_R_um, and Res_U_R_um
% save R(t) data in the selected time window t(tindstart:tindend)

Res_R_um_t_us      = zeros(length(t(tindstart:tindend)),4);
Res_R_um_t_us(:,1) = t(tindstart:tindend)*1e6; % all t data
Res_R_um_t_us(:,2) = R(tindstart:tindend)*1e6; % all R data
% Res_R_um_t_us(:,3) for selected t data and Res_R_um_t_us(:,4) for selected R data

% nr of curves that are to be visulized = nr. of all curves/Nsel
Nr_curves = floor(length(t(tindstart:tindend))/Nsel);

% selected p(r) curves
rsize = size(r_container);                         % the size of the array (T_KB,T_Gilmore)
Res_P_MPa_R_um = zeros(rsize(1,2)+10,2*Nr_curves); % rsize(1,2)=T_Gilmore, rsize(1,1)=T_KB,
                                                   % T_Gilmore +10 to include the point on the bubble wall and
                                                   % some points added at reshaping
                                                  
% selected u(r) curves                                                   
Res_U_R_um     = zeros(rsize(1,2)+10,2*Nr_curves); % rsize(1,2)=T_Gilmore, rsize(1,1)=T_KB,

% selected P,U,R values at bubble wall
Res_P0_U0_R0   = zeros(Nr_curves,3);            % save P U and R values at bubble wall
                                                   
% record the location of shock front peak for shock decay
Res_P_R_SW     = NaN(Nr_curves,2);              % save p(r),u(r) points at sw front                                               

%% time index transformation
% we have two time axles, one for bubble dynamics (adapative stepsize) and one
% for KB (fixed stepsize), we need to map the two time axis
ttind=ones(1,tindend);                       
for ii = tindstart:tindend
    ttind(ii) = floor((t(ii)-t(tindstart))/dtkb)+1;
end

%% shock wave reshaping and animation
%  simutaneous visulization of bubble dynamics and acoustic radiation,
%  subplot1 for bubble dynamics R(t) and subplot2 for acoustic radiation
hfig = figure();
set(hfig, 'units','normalized','outerposition',[0.05 0.05 0.9 0.9]);

h1 = subplot(1,2,1);    % subfig 1 R(t) curve, handle of axis
plot(t(tindstart:tindend)*1e6,R(tindstart:tindend)*1e6,'-b');
axis([min(t(tindstart:tindend)*1e6) max(t(tindstart:tindend)*1e6) 0 max(R(tindstart:tindend)*1e6*1.2)]);% x,y limits
xlabel('Time (\mus)');
ylabel('Bubble radius (\mum)');
set(gca,'fontsize',22);
hold on

h2=subplot(1,2,2);      % subfig 2 p(r) curve
loglog(Rtt(ttind(tindstart))*1e6,Ptt(ttind(tindstart))/1e6,'w');%'w' white, invisible
axis(xylim);            % set the limits of x and y
xlabel('Radius (\mum)','fontsize',fontsize)
ylabel('Pressure(MPa)','fontsize',fontsize)
set(gca,'fontsize',fontsize)   
hold on

% Set up the animation
Rstr  = [num2str(max(R)*1e6,2) 'um'];
tpstr = [num2str(tp,'%.2e') 's'];

if isARexp == 1
    Video_Name = ['Acoust_rad_at_expansion_' Rstr '_' tpstr '.avi'];
else
    Video_Name = ['Acoust_rad_at_collnreb_' Rstr '_' tpstr '.avi'];
end

% name of the output file in full path
toolboxPath = fileparts(fileparts(mfilename('fullpath'))); % Get the toolbox path
Video_Name  = fullfile(toolboxPath,'results',Video_Name);

writerObj = VideoWriter(Video_Name); 
writerObj.FrameRate = 5; % How many frames per second.
open(writerObj); 

for n_loop = tindstart:tindend      % go through the loops
    if rem(n_loop,Nsel) == 0        % select to display a curve by every "Nsel" curves, rem( a , b ) returns the remainder after division of a by b
        % shape the data: get rid of NaN and 0 values
        rplot = r_container(ttind(n_loop),:)*1e6;  %[1 x T_Gilmore] in um
        pplot = p_container(ttind(n_loop),:)/1e6;  %[1 x T_Gilmore] in MPa
        uplot = u_container(ttind(n_loop),:);      %[1 x T_Gilmore] m/s
        
        rplot(isnan(rplot)) = 0;            % we must get rid of the NaN's
        pplot(isnan(pplot)) = 0;            % we must get rid of the NaN's
        uplot(isnan(uplot)) = 0;            % we must get rid of the NaN's
        j = find(rplot);                    % and only find the valid values
        
        % fliping the array up side down to display progagating distance from
        % small to large, see the contour map of r(t_Gilmore,t_KB) in Figure 12 
        r_temp = flipud([rplot(j) Rtt(ttind(n_loop))*1e6]'); % [T_Gilmore x 1] in um transpose since x is in colume and after that flip up side down
        p_temp = flipud([pplot(j) Ptt(ttind(n_loop))/1e6]'); % [T_Gilmore x 1] in um transpose since y is in colume and after that flip up side down
        u_temp = flipud([uplot(j) Utt(ttind(n_loop))]');     % [T_Gilmore x 1]
        
        % save R P U values at bubble wall
        Res_P0_U0_R0(ceil((n_loop-tindstart+1)/Nsel),1) = Rtt(ttind(n_loop))*1e6; % save R value at bubble wall                                       
        Res_P0_U0_R0(ceil((n_loop-tindstart+1)/Nsel),2) = Ptt(ttind(n_loop))/1e6; % save P value at bubble wall  
        Res_P0_U0_R0(ceil((n_loop-tindstart+1)/Nsel),3) = Utt(ttind(n_loop));     % save U value at bubble wall  
        
        %% reshaping part I - bubble expansion phase
        if isARexp == 1
            % check the monotonic of curve
            N_monotonic = find(diff(r_temp)<0);
            if isempty(N_monotonic)  % the curve is monotonic and no shaping is needed, save y(x) directly
                r_reshaped = r_temp;
                p_reshaped = p_temp;
                u_reshaped = u_temp;
                loglog(r_temp, p_temp, 'b-','LineWidth',1.0);
                loglog(r_temp(1),p_temp(1),'ro');
                plot(h1,t(tindstart:tindend)*1e6,R(tindstart:tindend)*1e6,'-b','LineWidth',1.0);hold on
                plot(h1,t(n_loop)*1e6,R(n_loop)*1e6,'ro');hold on
                
                % save results
                Res_R_um_t_us(ceil((n_loop-tindstart+1)/Nsel),3) = t(n_loop)*1e6; % Save results, the red R(t) circles at selected time instants
                Res_R_um_t_us(ceil((n_loop-tindstart+1)/Nsel),4) = R(n_loop)*1e6; % Save results, the red R(t) circles at selected time instants
                Res_P_MPa_R_um(1:length(r_reshaped),2*ceil((n_loop-tindstart+1)/Nsel)-1) = r_reshaped;  % in um
                Res_P_MPa_R_um(1:length(p_reshaped),2*ceil((n_loop-tindstart+1)/Nsel))   = p_reshaped;  % in MPa
                Res_U_R_um(1:length(r_reshaped),2*ceil((n_loop-tindstart+1)/Nsel)-1) = r_reshaped;  % in um
                Res_U_R_um(1:length(p_reshaped),2*ceil((n_loop-tindstart+1)/Nsel))   = u_reshaped;  % in m/s
                
                Res_P_R_SW(ceil((n_loop-tindstart+1)/Nsel),1)  = NaN; % No sw front
                Res_P_R_SW(ceil((n_loop-tindstart+1)/Nsel),2)  = NaN;

            else % the curve is multi-valued and reshaping is required
               
                [r_reshaped,u_reshaped,p_reshaped,rout,pout] = shockFronttransfer(r_temp,u_temp,p_temp);
                
                % plot the reshaped p(r) curves
                loglog(r_reshaped, p_reshaped, 'b-','LineWidth',1.0);
                loglog(r_reshaped(1), p_reshaped(1), 'ro');
                % plot R(t) curve
                plot(h1,t(tindstart:tindend)*1e6,R(tindstart:tindend)*1e6,'-b','LineWidth',1.0);hold on
                plot(h1,t(n_loop)*1e6,R(n_loop)*1e6,'ro');hold on
                
                % save results 2020.05.02
                Res_R_um_t_us(ceil((n_loop-tindstart+1)/Nsel),3) = t(n_loop)*1e6; % Save results, the red R(t) circles at selected time instants
                Res_R_um_t_us(ceil((n_loop-tindstart+1)/Nsel),4) = R(n_loop)*1e6; % Save results, the red R(t) circles
                Res_P_MPa_R_um(1:length(r_reshaped),2*ceil((n_loop-tindstart+1)/Nsel)-1) = r_reshaped;  % in um
                Res_P_MPa_R_um(1:length(p_reshaped),2*ceil((n_loop-tindstart+1)/Nsel))   = p_reshaped;  % in MPa
                Res_U_R_um(1:length(r_reshaped),2*ceil((n_loop-tindstart+1)/Nsel)-1) = r_reshaped;  % in um
                Res_U_R_um(1:length(p_reshaped),2*ceil((n_loop-tindstart+1)/Nsel))   = u_reshaped;  % in m/s
                Res_P_R_SW(ceil((n_loop-tindstart+1)/Nsel),1)  = rout;   % record sw front position
                Res_P_R_SW(ceil((n_loop-tindstart+1)/Nsel),2)  = pout;
                
            end
        else
            %% reshaping part II - collapse and rebound phase
            if n_loop < tindRmin    % during collapse phase (before tRmin) curves are monotonic and no reshaping is needed - display and save data directly 
                r_reshaped = r_temp;
                p_reshaped = p_temp;
                u_reshaped = u_temp;
                
                loglog(h2,r_temp, p_temp,'r-','LineWidth',1.0);
                loglog(h2,r_temp(1),p_temp(1),'ro');
                plot(h1,t(tindstart:tindend)*1e6,R(tindstart:tindend)*1e6,'-b','LineWidth',1.0);hold on
                plot(h1,t(n_loop)*1e6,R(n_loop)*1e6,'ro');hold on
                
                % save results 
                Res_R_um_t_us(ceil((n_loop-tindstart+1)/Nsel),3) = t(n_loop)*1e6; % Save results, the red R(t) circles
                Res_R_um_t_us(ceil((n_loop-tindstart+1)/Nsel),4) = R(n_loop)*1e6; % Save results, the red R(t) circles
                Res_P_MPa_R_um(1:length(r_reshaped),2*ceil((n_loop-tindstart+1)/Nsel)-1) = r_reshaped;  % in um
                Res_P_MPa_R_um(1:length(p_reshaped),2*ceil((n_loop-tindstart+1)/Nsel))   = p_reshaped;  % in MPa
                Res_U_R_um(1:length(r_reshaped),2*ceil((n_loop-tindstart+1)/Nsel)-1) = r_reshaped;  % in um
                Res_U_R_um(1:length(p_reshaped),2*ceil((n_loop-tindstart+1)/Nsel))   = u_reshaped;  % in m/s
                Res_P_R_SW(ceil((n_loop-tindstart+1)/Nsel),1)  = NaN; 
                Res_P_R_SW(ceil((n_loop-tindstart+1)/Nsel),2)  = NaN;
                
            else  % during rebound phase curves are multi-valued and reshaping is required
                
                % check the monotonic 
                N_monotonic = find(diff(r_temp)<0);
                if isempty(N_monotonic)  % monotonic, no reshaping is needed
                    r_reshaped = r_temp;
                    p_reshaped = p_temp;
                    u_reshaped = u_temp;
                    loglog(h2,r_temp, p_temp,'b-','LineWidth',1.0);
                    loglog(h2,r_temp(1),p_temp(1),'bo');
                    plot(h1,t(tindstart:tindend)*1e6,R(tindstart:tindend)*1e6,'b-','LineWidth',1.0);hold on
                    plot(h1,t(n_loop)*1e6,R(n_loop)*1e6,'bo');hold on
                    
                    % save results 
                    Res_R_um_t_us(ceil((n_loop-tindstart+1)/Nsel),3) = t(n_loop)*1e6; % Save results, the red R(t) circles
                    Res_R_um_t_us(ceil((n_loop-tindstart+1)/Nsel),4) = R(n_loop)*1e6; % Save results, the red R(t) circles
                    Res_P_MPa_R_um(1:length(r_reshaped),2*ceil((n_loop-tindstart+1)/Nsel)-1) = r_reshaped;  % in um
                    Res_P_MPa_R_um(1:length(p_reshaped),2*ceil((n_loop-tindstart+1)/Nsel))   = p_reshaped;  % in MPa
                    Res_U_R_um(1:length(r_reshaped),2*ceil((n_loop-tindstart+1)/Nsel)-1) = r_reshaped;  % in um
                    Res_U_R_um(1:length(p_reshaped),2*ceil((n_loop-tindstart+1)/Nsel))   = u_reshaped;  % in m/s
                    Res_P_R_SW(ceil((n_loop-tindstart+1)/Nsel),1)  = NaN;
                    Res_P_R_SW(ceil((n_loop-tindstart+1)/Nsel),2)  = NaN;
                    
                else  % multi-valued and reshaping needed
                    [r_reshaped,u_reshaped,p_reshaped,rout,pout] = shockFronttransfer(r_temp,u_temp,p_temp);
                    
                    loglog(h2,r_reshaped, p_reshaped, 'b-','LineWidth',1.0);
                    loglog(h2,r_reshaped(1), p_reshaped(1), 'bo');
                    % plot rebound phase in blue
                    plot(h1,t(tindstart:tindend)*1e6,R(tindstart:tindend)*1e6,'b-','LineWidth',1.0);hold on
                    plot(h1,t(n_loop)*1e6,R(n_loop)*1e6,'bo');hold on
                    
                    % save results 
                    Res_R_um_t_us(ceil((n_loop-tindstart+1)/Nsel),3) = t(n_loop)*1e6; % Save results, the red R(t) circles
                    Res_R_um_t_us(ceil((n_loop-tindstart+1)/Nsel),4) = R(n_loop)*1e6; % Save results, the red R(t) circles
                    Res_P_MPa_R_um(1:length(r_reshaped),2*ceil((n_loop-tindstart+1)/Nsel)-1) = r_reshaped;  % in um
                    Res_P_MPa_R_um(1:length(p_reshaped),2*ceil((n_loop-tindstart+1)/Nsel))   = p_reshaped;  % in MPa
                    Res_U_R_um(1:length(r_reshaped),2*ceil((n_loop-tindstart+1)/Nsel)-1) = r_reshaped;  % in um
                    Res_U_R_um(1:length(p_reshaped),2*ceil((n_loop-tindstart+1)/Nsel))   = u_reshaped;  % in m/s
                    Res_P_R_SW(ceil((n_loop-tindstart+1)/Nsel),1)  = rout;
                    Res_P_R_SW(ceil((n_loop-tindstart+1)/Nsel),2)  = pout;
                end
                
            end
            
            
        end
               
        frame = getframe(gcf); 
        writeVideo(writerObj, frame);
        %pause(0.001);
    end
end

hold off
close (writerObj); % save the movie
disp('----                   Animation complete                   ----')
%% replacing zeros by NaN for the sake of plotting in software like OriginLab
[row,col] = find(Res_R_um_t_us==0); 
for i=1:length(row)
    Res_R_um_t_us(row(i),col(i)) =NaN;
end

[row,col] = find(Res_P_MPa_R_um==0);
for i=1:length(row)
    Res_P_MPa_R_um(row(i),col(i)) = NaN;
end


%% plot Nr_curves of u(r) and p(r)

fontsize = 14;
if isExtrafig == 1
figure() % u(r) curves
for i=1:Nr_curves
semilogx(Res_P0_U0_R0(:,1),Res_P0_U0_R0(:,3),'ro',Res_P0_U0_R0(:,1),Res_P0_U0_R0(:,3),'b--',Res_U_R_um(:,2*i-1),Res_U_R_um(:,2*i),'b-','LineWidth',1.0)
hold on
end
hold off
xlabel('Radius (\mum)','fontsize',fontsize)
ylabel('Velocity(m/s)','fontsize',fontsize)
set(gca,'fontsize',fontsize)

figure() % p(r) curves
for i=1:Nr_curves
loglog(Res_P0_U0_R0(:,1),Res_P0_U0_R0(:,2),'ro',Res_P0_U0_R0(:,1),Res_P0_U0_R0(:,2),'b--',Res_P_R_SW(:,1),Res_P_R_SW(:,2),'k-.',Res_P_MPa_R_um(:,2*i-1),Res_P_MPa_R_um(:,2*i),'b-','LineWidth',1.0)
hold on
end
hold off
xlabel('Radius (\mum)','fontsize',fontsize)
ylabel('Pressure(MPa)','fontsize',fontsize)
set(gca,'fontsize',fontsize)
end

%% plot the finely selected curves for publication
% The purpose of this part is to select and store curves directly in
% Matlab such that tedious work in OriginLab can be avoided and Figures as Figs.9 and 10 in Liang2022JFM can be created directly.

Nsel_t        = Res_R_um_t_us(:,3)*1e3; % Nsel - number of selected t curves in unit ns  
Nsel_t(isnan(Nsel_t)) = 0;              % get rid of NaNs
j             = find(Nsel_t);           % find valid values
Nsel_t        = Nsel_t(j);              % resize array
sel_ind_array = interp1(Nsel_t, 1:length(Nsel_t), time_q,'nearest'); % time index of the querried time in array Nsel_t
assert (max(sel_ind_array) <= Nr_curves && min(sel_ind_array)>0);    % make sure the index is within the range of nr_datapair

sel_ind_array(isnan(sel_ind_array)) = 0;% get rid of NaNs
j             = find(sel_ind_array);    % find valid values
sel_ind_array = sel_ind_array(j);       % resize array
% sel_ind_array = [1 3 5 9 12 15 18 21 24 32 45 71 Nr_curves]; % use index directly
nr = length(sel_ind_array); % nr of selected curves

% preallocation of arrays for saving results of selected p(r) and u(r) curves
Res_Sel_R_P      = zeros(length(Res_P_MPa_R_um),2*nr); % R in um and P in MPa 
Res_Sel_R_U      = zeros(length(Res_U_R_um),2*nr);     % R in um and U in m/s
Res_Sel_R0_P0_U0 = zeros(nr,3);                        % R0, P0 and U0 
Res_Sel_P_R_SW   = zeros(nr,2);                        % R_sw, P_sw
Res_Sel_t_us_R_um= zeros(nr,2);                        % t, R0

for i = 1:nr
    index = sel_ind_array(i);
    Res_Sel_R_P(:,2*i-1) = Res_P_MPa_R_um(:,2*index-1); % save slected r
    Res_Sel_R_P(:,2*i)   = Res_P_MPa_R_um(:,2*index);   % save slected p
    
    Res_Sel_R_U(:,2*i-1) = Res_U_R_um(:,2*index-1);     % save slected r
    Res_Sel_R_U(:,2*i)   = Res_U_R_um(:,2*index);       % save slected u
    
    Res_Sel_R0_P0_U0(i,1) = Res_P0_U0_R0(index,1);      % save slected R0
    Res_Sel_R0_P0_U0(i,2) = Res_P0_U0_R0(index,2);      % slected P0
    Res_Sel_R0_P0_U0(i,3) = Res_P0_U0_R0(index,3);      % slected U0
    
    Res_Sel_P_R_SW (i,1)  = Res_P_R_SW(index,1);        % slected r_sw
    Res_Sel_P_R_SW (i,2)  = Res_P_R_SW(index,2);        % slected p_sw
    
    Res_Sel_t_us_R_um(i,1)= Res_R_um_t_us(index,3);     % slected t in us  
    Res_Sel_t_us_R_um(i,2)= Res_R_um_t_us(index,4);     % slected R0 in um
end

figure() % selected 10 p(r) curves
loglog(Res_P0_U0_R0(:,1),Res_P0_U0_R0(:,2),'b--',...        % all points at bubble wall
       Res_Sel_R0_P0_U0(:,1),Res_Sel_R0_P0_U0(:,2),'ro',... % selected points at bubble wall
       Res_P_R_SW(:,1),Res_P_R_SW(:,2),'b-.',...            % all points at shock front
       Res_Sel_P_R_SW(:,1),Res_Sel_P_R_SW(:,2),'bs','LineWidth',1.0) % selected points at shock front
hold on
for i=1:nr
    loglog(Res_Sel_R_P(:,2*i-1),Res_Sel_R_P(:,2*i),'k','LineWidth',1.0)
     txt1 = ['\leftarrow ' num2str(Res_Sel_t_us_R_um(i,1)*1e3,'%.2f') ' ns'];
     text(Res_Sel_R0_P0_U0(i,1),Res_Sel_R0_P0_U0(i,2),txt1,'HorizontalAlignment','left', 'FontName', 'Times New Roman');
end
hold off
xlabel('Radius (μm)','fontsize',fontsize)
ylabel('Pressure(MPa)','fontsize',fontsize)
set(gca,'fontsize',fontsize)

figure() % selected 10 u(r) curves
semilogx(Res_P0_U0_R0(:,1),Res_P0_U0_R0(:,3),'b--',Res_Sel_R0_P0_U0(:,1),Res_Sel_R0_P0_U0(:,3),'ro','LineWidth',1.0)
hold on
for i=1:nr
    semilogx(Res_Sel_R_U(:,2*i-1),Res_Sel_R_U(:,2*i),'k','LineWidth',1.0)
     txt1 = ['\leftarrow ' num2str(Res_Sel_t_us_R_um(i,1)*1e3,'%.1f') ' ns'];
     text(Res_Sel_R0_P0_U0(i,1),Res_Sel_R0_P0_U0(i,3),txt1,'HorizontalAlignment','left', 'FontName', 'Times New Roman');
end
hold off
xlabel('Radius (μm)','fontsize',fontsize)
ylabel('Velocity(m/s)','fontsize',fontsize)
set(gca,'fontsize',fontsize)

%% Plot shock front radus-time and shock velocity-time curves
B = [Res_Sel_t_us_R_um(:,1)*1e3,Res_Sel_P_R_SW(:,1)];
[~, ia, ~] = unique(B(:, 1), 'first');
SW_radius_time_curve = B(ia, :);

% figure() 
% plot(SW_radius_time_curve(:,1),SW_radius_time_curve(:,2),'b','LineWidth',1)
% xlabel('Time (ns)','fontsize',fontsize)
% ylabel('Shock Wave Radius (μm)','fontsize',fontsize)
% set(gca,'fontsize',fontsize)

SW_velocity_time_curve = zeros(length(SW_radius_time_curve(:,1))-1,2);
for i = 1:1:length(SW_radius_time_curve(:,1))-1
    tt = (SW_radius_time_curve(i+1,1)-SW_radius_time_curve(i,1))/2+SW_radius_time_curve(i,1);
    uu = (SW_radius_time_curve(i+1,2)-SW_radius_time_curve(i,2))/(SW_radius_time_curve(i+1,1)-SW_radius_time_curve(i,1))*1000;
    SW_velocity_time_curve(i,:) = [tt,uu];
end

% figure() 
% plot(SW_velocity_time_curve(:,1),SW_velocity_time_curve(:,2),'b','LineWidth',1)
% xlabel('Time (ns)','fontsize',fontsize)
% ylabel('Shock Wave Velocity (μm)','fontsize',fontsize)
% set(gca,'fontsize',fontsize)

%% Save results for outputs
varargout{1} = Res_Sel_R_P;
varargout{2} = Res_Sel_R_U;
varargout{3} = Res_Sel_R0_P0_U0;
varargout{4} = Res_Sel_P_R_SW; % P(r) curve for selected shock fronts
varargout{5} = Res_Sel_t_us_R_um;
varargout{6} = SW_radius_time_curve;
varargout{7} = SW_velocity_time_curve;
varargout{8} = Res_P_R_SW; % P(r) curve for all calculated shock fronts
varargout{9} = Res_P0_U0_R0; %

end