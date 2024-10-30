close all
clearvars

%% Show test header (useful when run as part of test suite)
disp('-----------------------------------------------------------------')
disp('-----            LIBDAR TOOLBOX TEST:                       -----')
disp('-----      Bubble dynamics and acoustic radiation example   -----')
disp('-----   Ref: Wen et al. Ultrason Sonochem 95, 106391 (2023) -----')
disp('Acknowledge:Thank Prof. Z.F. Yao for sharing the exp. R(t) data.')
disp(' ')

%% Add path to necessary functions
toolboxPath=fileparts(fileparts(mfilename('fullpath'))); %Get the toolbox path

%Add core, misc and expData path
addpath(fullfile(toolboxPath,'core_functions'),fullfile(toolboxPath,'misc'),fullfile(toolboxPath,'expData'));

%% explicit input parameters
tp   = 5e-9;      % pulse duration, 5 ns in the Ref. Wen2023UltrasonSonochem
T1   = 222.31e-6; % measured 1st osc. time, 222.31 us, see Table 2 in the Ref. Wen2023 
T2   = 88.22e-6;  % measured 2nd osc. time, 88.22 us, see Table 2 in the Ref. Wen2023
R0   = 19e-6;     % initial bubble radius, taken as exp. determined plasma radius, see Fig. 7 in the Ref. Wen2023
Rnbd = 439.7e-6;  % equilibrium bubble radius after bd, tuning parameter such that simulated Tosc1 fits with measured value
                  % Note: Rnbd value is slightly larger than that in the Ref. Wen2023 because here we use the general "jump-start" condition, 
                  % while in the Ref. Wen2023 2nd-order jump-start condition under inertial-confinement was used.
Rnc1 = 187e-6;    % equilibrium bubble radius at 1st bubble collapse, tuning parameter such that simulated Tosc2 fits with measured value
Rnc2 = 107e-6;    % equilibrium bubble radius at 2nd bubble collapse, tuning parameter such that simulated Tosc3 fits with measured value

disp('Simulating...')
tic

%% load experimental data
bubbledata = xlsread('expData_1.2mm_5ns.xlsx');
time_exp   = bubbledata(:,1); % in us
radius_exp = bubbledata(:,2); % in um

%% bubble dynamics


[filename,t,R,P,U,tind_sw,tindRmin] = bubble_dynamics_and_energy_partitioning (T1,T2,R0,Rnbd,Rnc1,Rnc2,tp,'GS','adaptive','MaxStepSize',2e-9,'EP',0,'JSC','general','tRmax3',336e-6,'Rnc3',82.5e-6,'isIsotherm',0,'tend',400e-6);

% exp. v.s. cal. curves
figure();
plot(time_exp, radius_exp, 'ro',t*1e6,R*1e6,'b-','LineWidth',1.5);
xlabel('Time (\mus)')
ylabel('Radius (\mum)')
legend('exp.','theo')
title('R(t) exp. v.s. theo.');
set(gca,'fontsize',14)

%% acoustic radiation
isARexp = 1; % is acoustic radiation at bubble expansion? 1 for bubble expansion and 0 for collapse and rebound

if isARexp == 1 % acoustic radiation during bubble expansion     

    [tindstart,tindend,dtkb,r_container,u_container,p_container,Rtt,Utt,Ptt] = KB_acoustic_radiation(t,R,P,U,tind_sw,'isARexp',isARexp,'winSel','manuel','dtkb',2e-12,'t_sw0',t(2),'t_swd',500e-9);

    % shock front before and after reshaping at queried times, reproducing Fig. S1 in Liang2022JFM
    time_q = [10,15,20]; % querried time in unit ns, 
    [varargout] = SelcurvReshaping(t,tindstart,tindend,dtkb,r_container,u_container,p_container,Rtt,Utt,Ptt,'time_q',time_q);

    % acoustic radiation, animation, and publication quality plots, reproducing Figs. 9 and 10 in Liang2022JFM
    time_q_pub = [5,6,7,8,9,10,12,15,20,30,50,90,160,280,400];% querried time in unit ns
    acoustRadwReshaping(t,R,tp,tindstart,tindend,tindRmin,dtkb,r_container,u_container,p_container,Rtt,Utt,Ptt,'time_q',time_q_pub,'isARexp',isARexp,'xylim',[10 1000 0.1 1e4],'Nsel',8,'isExtrafig',0);

else  % acoustic radiation during bubble collapse and rebound   
    tcoll = 222.173e3; % time at bubble collapse in unit ns, see Fig. 9 in the Ref. Wen2023
    [tindstart,tindend,dtkb,r_container,u_container,p_container,Rtt,Utt,Ptt] = KB_acoustic_radiation(t,R,P,U,tind_sw,'isARexp',isARexp,'winSel','manuel','dtkb',2e-12,'t_sw0',tcoll*1e-9-150e-9,'t_swd',tcoll*1e-9+220e-9);

    % shock front before and after reshaping at queried times, reproducing Fig. S1 in Liang2022JFM
    time_q = tcoll+[10,15,20]; % querried time in unit ns
    [varargout] = SelcurvReshaping(t,tindstart,tindend,dtkb,r_container,u_container,p_container,Rtt,Utt,Ptt,'time_q',time_q);

    % acoustic radiation, animation, and publication quality plots, reproducing Figs. 9 and 10 in Liang2022JFM
    time_q_pub = tcoll+[-150,-100,-60,-30,-15,-8,-4,-2,-1,0,1,2,4,8,15,30,60,100,150,220]; % querried time in unit ns
    acoustRadwReshaping(t,R,tp,tindstart,tindend,tindRmin,dtkb,r_container,u_container,p_container,Rtt,Utt,Ptt,'time_q',time_q_pub,'isARexp',isARexp,'xylim',[10 1000 1 6000],'Nsel',8,'isExtrafig',0);

end

toc

%% Add blank line (nicer formatting for test text output).
disp('-----                  End of simulation                -------- ')