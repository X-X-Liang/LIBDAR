close all
clearvars

%% Show test header (useful when run as part of test suite)
disp('-----------------------------------------------------------------')
disp('-----            LIBDAR TOOLBOX TEST:                       -----')
disp('-----    Bubble dynamics and energy partitioning example    -----')
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
Rnbd = 13.718e-6; % equilibrium bubble radius after breakdown, tuning parameter such that simulated Tosc1 fits with measured value, see Table 1 in the Ref. Liang2022JFM
Rnc1 = 3.615e-6;  % equilibrium bubble radius at 1st bubble collapse, tuning parameter such that simulated Tosc2 fits with measured value, see Table 1 in the Ref. Liang2022JFM
Rnc2 = 2.415e-6;  % equilibrium bubble radius at 2nd bubble collapse, tuning parameter such that simulated Tosc3 fits with measured value, see Table 1 in the Ref. Liang2022JFM

%% optional parameters for varargin with name-value combinations
% 'JSC' - jump-start conditions with options 'no','first-order','second-order','general'(default)
% 'GS'  - Gilmore Solver with options, 'adaptive'(default), 'maxstepsize', 'fixedstepsize'
% 'EP'  - calculating energy partitioning, 1 for true (default), 0 for false
% 'tend'- simulation ends at tend, default value tend = 1.5*(T1+T2)
% 'tRmax3'      - time instant at Rmax3
% 'Twall'       - temperature at the bubble wall
% 'StepSize'    - tuning step size when GS = 'fixedstepsize'is selected
% 'MaxStepSize' - tuning max. step size when GS = 'maxstepsize'is selected

%% call the function
disp('Simulating...')
tic
% Reproducing Figs. 2, 7, 8, 12 and Table 2 in the Ref: Liang2022JFM using the following command
% Results on energy partitioning is stored in the excel file 'filename.xlsx'

[filename,t,R,pg,U,tind_sw] = bubble_dynamics_and_energy_partitioning (T1,T2,R0,Rnbd,Rnc1,Rnc2,tp,'JSC','general','GS','adaptive','EP',1);

% Reproducing Fig. 11 in the Ref: Liang2022JFM using the following command
% [filename, ~] = bubble_dynamics_and_energy_partitioning (T1,T2,R0,Rnbd,Rnc1,Rnc2,tp,'JSC','general','GS','adaptive','EP',0,'tend',90e-6,'tRmax3',8.5203e-6,'Twall',110+273);
toc

%% Add blank line (nicer formatting for test text output).
disp('-----                  End of simulation                 -------- ')