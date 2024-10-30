close all
clearvars

%% Show test header (useful when run as part of test suite)
disp('-----------------------------------------------------------------')
disp('-----            LIBDAR TOOLBOX TEST:                       -----')
disp('-----      Bubble dynamics and acoustic radiation example   -----')
disp('-----    Ref: Liang et al. J. Fluid. Mech. 940, A5 (2022)   -----')
disp('-----------------------------------------------------------------')
disp(' ')

%% Add path to necessary functions
toolboxPath=fileparts(fileparts(mfilename('fullpath'))); %Get the toolbox path

%Add core and misc path
addpath(fullfile(toolboxPath,'core_functions'),fullfile(toolboxPath,'misc'),fullfile(toolboxPath,'results'));

%% explicit input parameters
tp   = 265e-15;   % pulse duration, 265e-15 in Liang2022JFM 
T1   = 6.488e-6;  % measured 1st osc. time, 6.488e-6, see Fig. 6 and Table 1 in the Ref. Liang2022JFM 
T2   = 1.5616e-6; % measured 2nd osc. time, 1.562e-6, see Fig. 6 and Table 1 in the Ref. Liang2022JFM
R0   = 1.33e-6;   % initial bubble radius, taken as exp. determined plasma radius, see Fig. 4 in the Ref. Liang2022JFM
Rnbd = 13.718e-6; % equilibrium bubble radius after breakdown, tuning parameter such that simulated Tosc1 fits with measured value
Rnc1 = 3.615e-6;  % equilibrium bubble radius at 1st bubble collapse, tuning parameter such that simulated Tosc2 fits with measured value
Rnc2 = 2.415e-6;  % equilibrium bubble radius at 2nd bubble collapse, tuning parameter such that simulated Tosc3 fits with measured value

%% optional parameters for bubble dynamics and energy partitioning
% 'JSC' - jump-start conditions with options 'no','first-order','second-order','general'(default)
% 'GS'  - Gilmore Solver with options, 'adaptive'(default), 'maxstepsize', 'fixedstepsize'. 
%         'adaptive' for bubble dynamics and 'maxstepsize' for acoustic radiation
% 'EP'  - calculating energy partitioning, 1 for true (default), 0 for false
% 'tend'- end time for bubble dynamics simulation, default value tend = 1.5*(T1+T2)
% 'tRmax3'      - time instant at Rmax3
% 'Twall'       - temperature at the bubble wall
% 'MaxStepSize' - tuning max. step size when GS = 'maxstepsize'is selected

%% optional parameters for acoustic emission
% 't_sw0'   - start time for shock wave simulation
% 't_swd'   - end time for shock wave simulation 
% 'dtkb'    - fixed time stepsize for kickwood-bethe calculation, which depends on Rmax 
%             % for Rmax ~ 1µm, try dtkb = 400fs
%             % for 10 µm < Rmax < 100 µm, try dtkb = 2 ps, 
%             % for Rmax ~ 500 µm, try dtkb = 20 ps 
% 'isARexp' - choice of shock wave phase. Options: 
%             1 for bubble expansion and 0 for bubble collapse and rebound
% 'winSel'  - choice of selecting simulation time window. Options:
%             'cursor' - using the mouse cursor with visual interaction
%             'manual' - manually give time instants t_sw0 and t_swd

%% optional parameters for shock wave reshaping and visulization
% 'time_q' - queried time instants for visulization, in unit nanosecond

disp('Simulating...')
tic
%% bubble dynamics

[filename,t,R,P,U,tind_sw,tindRmin] = bubble_dynamics_and_energy_partitioning (T1,T2,R0,Rnbd,Rnc1,Rnc2,tp,'GS','maxstepsize','MaxStepSize',0.3e-9,'EP',0);

%% acoustic radiation
isARexp = 0; % is acoustic radiation at bubble expansion? 1 for bubble expansion and 0 for collapse and rebound

if isARexp == 1 % acoustic radiation during bubble expansion     

    [tindstart,tindend,dtkb,r_container,u_container,p_container,Rtt,Utt,Ptt] = KB_acoustic_radiation(t,R,P,U,tind_sw,'isARexp',isARexp,'winSel','manuel','dtkb',4e-12,'t_sw0',t(2),'t_swd',135e-9);
     
    % shock front before and after reshaping at queried times, reproducing Fig. S1 in Liang2022JFM
    time_q = [0.33,1.64,8.72]; % querried time in unit ns, 
    [varargout] = SelcurvReshaping(t,tindstart,tindend,dtkb,r_container,u_container,p_container,Rtt,Utt,Ptt,'time_q',time_q);

    % acoustic radiation, animation, and publication quality plots, reproducing Figs. 9 and 10 in Liang2022JFM
    time_q_pub = [0.05,0.7,2,4.7,9.7,19.7,37.7,70.7,130];% querried time in unit ns
    acoustRadwReshaping(t,R,tp,tindstart,tindend,tindRmin,dtkb,r_container,u_container,p_container,Rtt,Utt,Ptt,'time_q',time_q_pub,'isARexp',isARexp,'xylim',[1 300 0.1 3e3],'Nsel',4,'isExtrafig',0);

else  % acoustic radiation during bubble collapse and rebound   

    [tindstart,tindend,dtkb,r_container,u_container,p_container,Rtt,Utt,Ptt] = KB_acoustic_radiation(t,R,P,U,tind_sw,'isARexp',isARexp,'winSel','manuel','dtkb',4e-12,'t_sw0',6480e-9,'t_swd',6515e-9);
    
    % shock front before and after reshaping at queried times, reproducing Fig. S1 in Liang2022JFM
    tcoll = 6488; % time at bubble collapse in unit ns, see Fig. 10 in Liang2022JFM
    time_q = tcoll+[0.054,0.155,0.397]; % querried time in unit ns
    [varargout] = SelcurvReshaping(t,tindstart,tindend,dtkb,r_container,u_container,p_container,Rtt,Utt,Ptt,'time_q',time_q);
    
    % acoustic radiation, animation, and publication quality plots, reproducing Figs. 9 and 10 in Liang2022JFM
    time_q_pub = tcoll+[-1.15,-0.57,-0.27,-0.13,-0.07,-0.03,0,+0.09,+0.3,+0.76,+1.45,+3.2,+7.03,+14.03]; % querried time in unit ns
    acoustRadwReshaping(t,R,tp,tindstart,tindend,tindRmin,dtkb,r_container,u_container,p_container,Rtt,Utt,Ptt,'time_q',time_q_pub,'isARexp',isARexp,'xylim',[0.2 70 0.1 3e4],'Nsel',4,'isExtrafig',0);

end

toc

%% Add blank line (nicer formatting for test text output).
disp('-----                  End of simulation                 ------- ')
disp(' ')