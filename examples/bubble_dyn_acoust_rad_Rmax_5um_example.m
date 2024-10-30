close all
clearvars

%% Show test header (useful when run as part of test suite)
disp('-----------------------------------------------------------------')
disp('-----            LIBDAR TOOLBOX TEST:                       -----')
disp('-----      Bubble dynamics and acoustic radiation example   -----')
disp('-----       Ref: Liang&Vogel Preprint (2024)                -----')
disp('-----------------------------------------------------------------')
disp(' ')

%% Add path to necessary functions
toolboxPath=fileparts(fileparts(mfilename('fullpath'))); %Get the toolbox path

%Add core and misc path
addpath(fullfile(toolboxPath,'core_functions'),fullfile(toolboxPath,'misc'),fullfile(toolboxPath,'results'));

%% explicit input parameters
tp   = 100e-12;  % pulse duration, 100 fs, 100 ps, 500 ps, and 5 ns
T1   = 0.9e-6;   % 1st osc. 0.9136 us for 100 fs; and 0.9 us for 100 ps,0.87 us for 500 ps and 0.74 us for 5 ns
T2   = 0.1551e-6;% 2nd osc. 
R0   = 210e-9;   % initial bubble radius
Rnbd = 10.4*R0;  % equilibrium bubble radius after breakdown, 
Rnc1 = Rnbd/4;   % equilibrium bubble radius at 1st bubble collapse, 
                 % Note: for studying bubble expansion, the value of Rnc is not important. Here we use a ratio value close to that in the Ref. Liang2022 JFM. 
                 % In reality, Rnc depends on bubble size as well as plasma energy density.
Rnc2 = Rnbd/4;   % equilibrium bubble radius at 2nd bubble collapse,

disp('Simulating...')
tic
%% bubble dynamics

[filename,t,R,P,U,tind_sw,tindRmin] = bubble_dynamics_and_energy_partitioning (T1,T2,R0,Rnbd,Rnc1,Rnc2,tp,'GS','maxstepsize','MaxStepSize',0.2e-9,'EP',0);

%% acoustic radiation
isARexp = 1; % is acoustic radiation at bubble expansion? 1 for bubble expansion and 0 for collapse and rebound

if isARexp == 1 % acoustic radiation during bubble expansion     

    [tindstart,tindend,dtkb,r_container,u_container,p_container,Rtt,Utt,Ptt] = KB_acoustic_radiation(t,R,P,U,tind_sw,'isARexp',isARexp,'winSel','manuel','dtkb',4e-12,'t_sw0',t(2),'t_swd',16e-9);
     
    % acoustic radiation, animation, and publication quality plots, reproducing Figs.12&13 in Liang&Vogel2024 preprint
    time_q_pub = [0.001,0.036,0.14,0.21,0.25,0.4,1,1.8,4.4,9.6,11,15]; % querried time in unit ns
    acoustRadwReshaping(t,R,tp,tindstart,tindend,tindRmin,dtkb,r_container,u_container,p_container,Rtt,Utt,Ptt,'time_q',time_q_pub,'isARexp',isARexp,'xylim',[0.1 40 0.1 3e3],'Nsel',4,'isExtrafig',0);

else  % acoustic radiation during bubble collapse and rebound   

    % since we study the expansion phase only, the collapse and rebound is
    % empty here because Tosc2 is not measured and the info on Rnc1 is not accurate

end

toc

%% Add blank line (nicer formatting for test text output).
disp('-----                  End of simulation                 -------- ')