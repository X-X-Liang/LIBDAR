close all
clearvars

%% Show test header (useful when run as part of test suite)
disp('-----------------------------------------------------------------')
disp('-----            LIBDAR TOOLBOX TEST:                       -----')
disp('-----      Bubble dynamics and acoustic radiation example   -----')
disp('-----   Ref: Interferometry data 355nm, 560ps, Rmax=15.3um  -----')
disp(' ')

%% Add path to necessary functions
toolboxPath=fileparts(fileparts(mfilename('fullpath'))); %Get the toolbox path

%Add core, misc and expData path
addpath(fullfile(toolboxPath,'core_functions'),fullfile(toolboxPath,'misc'),fullfile(toolboxPath,'expData'));

%% explicit input parameters
tp   = 0.56e-9;   % pulse duration, 0.56 ns 
T1   = 2691e-9; % measured 1st osc. time, 222.31 us, see Table 2 in the Ref. Wen2023 
T2   = 907.9e-9;  % measured 2nd osc. time, 88.22 us, see Table 2 in the Ref. Wen2023
R0   = 0.433e-6;     % initial bubble radius, taken as exp. determined plasma radius, see Fig. 7 in the Ref. Wen2023
Rnbd = 13.657*R0;  % equilibrium bubble radius after bd, tuning parameter such that simulated Tosc1 fits with measured value
                  % Note: Rnbd value is slightly larger than that in the Ref. Wen2023 because here we use the general "jump-start" condition, 
                  % while in the Ref. Wen2023 2nd-order jump-start condition under inertial-confinement was used.
Rnc1 = 5.165*R0;    % equilibrium bubble radius at 1st bubble collapse, tuning parameter such that simulated Tosc2 fits with measured value
Rnc2 = 1.87*R0;    % equilibrium bubble radius at 2nd bubble collapse, tuning parameter such that simulated Tosc3 fits with measured value

disp('Simulating...')
tic

%% load experimental data
bubbledata = xlsread('expData_15.3um_0.56ns_adjusted.xlsx');
time_exp   = bubbledata(:,1); % in us
radius_exp = bubbledata(:,2); % in um

[Rmax,i_Rmax] = max(radius_exp);
avg_t_exp = time_exp(1:i_Rmax);
avg_R_exp = avg_t_exp;
Delta_t=avg_t_exp;
Delta_R=avg_t_exp;
U_exp=avg_t_exp;

for i=1:(i_Rmax-1)
    avg_t_exp(i) = (time_exp(i)+time_exp(i+1))/2;
    avg_R_exp(i) = (avg_R_exp(i)+avg_R_exp(i+1))/2;
    Delta_t(i)=time_exp(i+1)-time_exp(i);
    Delta_R(i)=avg_R_exp(i+1)-avg_R_exp(i);
    U_exp(i) = Delta_R(i)/Delta_t(i);
end


%% bubble dynamics

% [filename,t,R,P,U,tind_sw,tindRmin] = bubble_dynamics_and_energy_partitioning (T1,T2,R0,Rnbd,Rnc1,Rnc2,tp,'GS','adaptive','MaxStepSize',2e-9,'EP',0,'JSC','general','tRmax3',336e-6,'Rnc3',82.5e-6,'isIsotherm',0,'tend',400e-6);
[filename,t,R,P,U,tind_sw,tindRmin] = bubble_dynamics_and_energy_partitioning (T1,T2,R0,Rnbd,Rnc1,Rnc2,tp,'GS','adaptive','EP',0,'JSC','general');

% exp. v.s. cal. curves
figure();
plot(time_exp, radius_exp, 'ro',t*1e6,R*1e6,'b-','LineWidth',1.5);
xlabel('Time (\mus)')
ylabel('Radius (\mum)')
legend('exp.','theo')
title('R(t) exp. v.s. theo.');
set(gca,'fontsize',14)

figure();
plot(avg_t_exp, U_exp, 'ro',t*1e6,U,'b-','LineWidth',1.5);
xlabel('Time (\mus)')
ylabel('Velocity (\mum)')
legend('exp.','theo')
title('U(t) expantion');
set(gca,'fontsize',14)

toc

%% Add blank line (nicer formatting for test text output).
disp('-----                  End of simulation                -------- ')