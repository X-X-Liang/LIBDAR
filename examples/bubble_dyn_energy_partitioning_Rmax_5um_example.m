close all
clearvars

%% Show test header (useful when run as part of test suite)
disp('-----------------------------------------------------------------')
disp('-----            LIBDAR TOOLBOX TEST:                       -----')
disp('-----    Bubble dynamics and energy partitioning example    -----')
disp('-----       Ref: Liang&Vogel Preprint (2024)                -----')
disp('-----------------------------------------------------------------')
disp(' ')

%% Add path to necessary functions
toolboxPath=fileparts(fileparts(mfilename('fullpath'))); %Get the toolbox path

%Add core and misc path
addpath(fullfile(toolboxPath,'core_functions'),fullfile(toolboxPath,'misc'),fullfile(toolboxPath,'results'));

%% explicit input parameters
tp   = 5e-9;  % pulse duration, 100 fs, 100 ps, 500 ps, and 5 ns
T1   = 0.74e-6;   % 1st osc. 0.9136 us for 100 fs; and 0.9 us for 100 ps,0.87 us for 500 ps and 0.74 us for 5 ns
T2   = 0.155e-6; % 2nd osc. 
R0   = 210e-9;   % initial bubble radius
Rnbd = 10.4*R0;  % equilibrium bubble radius after breakdown, 
Rnc1 = Rnbd/4;   % equilibrium bubble radius at 1st bubble collapse, 
                 % Note: for studying bubble expansion, the value of Rnc is not important. Here we use a ratio value close to that in the Ref. Liang2022 JFM. 
                 % In reality, Rnc depends on bubble size as well as plasma energy density.
Rnc2 = Rnbd/4;   % equilibrium bubble radius at 2nd bubble collapse,

%% call the function
disp('Simulating...')
tic
% Reproducing Fig. 9 in the Ref: Liang&Vogel Preprint (2024) using the following command
% Results on energy partitioning is stored in the excel file 'filename.xlsx'
% Tuning R0 and pulse duration, and reading energy partitioning data from filename.xlsx will create Figs.8, 10 and 11 in the Ref: Liang&Vogel Preprint (2024). 

[filename,t,R,pg,U,tind_sw] = bubble_dynamics_and_energy_partitioning (T1,T2,R0,Rnbd,Rnc1,Rnc2,tp,'JSC','general','GS','adaptive','EP',1);


toc

%% Add blank line (nicer formatting for test text output).
disp('-----                  End of simulation                 -------- ')