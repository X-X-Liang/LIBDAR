close all
clearvars

%% Show test header (useful when run as part of test suite)
disp('-----------------------------------------------------------------')
disp('-----            LIBDAR TOOLBOX TEST:                       -----')
disp('-----  Inertial confinement (Rmax,tau_L) example            -----')
disp('----Ref: Liang&Vogel,https://doi.org/10.48550/arXiv.2501.13749---')
disp('-----------------------------------------------------------------')
disp(' ')

%% Add path to necessary functions
toolboxPath=fileparts(fileparts(mfilename('fullpath'))); %Get the toolbox path

%Add core and misc path
addpath(fullfile(toolboxPath,'core_functions'),fullfile(toolboxPath,'misc'));

%% Input parameters
Rnbd_R0_ratio = 10.4;           % To reproduce the results in the Ref., one can set here 5, 10.4 and 15
Inertial_threshold = 2^(1/3);   % Definition of inertial confinement, see Eq. (1) in the Ref.

%% Define parameter space (tau_L,R0)
NoP = 21; % number of points, 51 points for fine grids - computationally slow; 21 points for coarse grids- computationally fast
tau_Array = logspace(-13,-8,NoP);  % array of pulse duration, generates NoP points between 100 fs 10^-13 and 10 ns 10^-8.
L_tau     = length(tau_Array);

R0_Array  = logspace(-8.2,-4,NoP); % array of plasma size, generates NoP points between ~ 10 nm and 0.1 mm.
                                  % Note: To reproduce Fig. 6 in the Ref., set R0_Array  = 210e-9; 
                                  % Note: To reproduce Fig. 3(a) and 4 in the Ref., set
                                  % R0_Array=logspace(-7.9:-4,51) for Rnbd_R0_ratio=5; 
                                  % R0_Array=logspace(-8.2:-4,51) for Rnbd_R0_ratio=10.4;
                                  % R0_Array=logspace(-8.4:-4,51) for Rnbd_R0_ratio=15; 
L_R0      = length(R0_Array);

% Initialization of arrays and matrix to store results
tau_th_Array  = zeros(L_R0,1);
Rmax_th_Array = zeros(L_R0,1);
Tosc_th_Array = zeros(L_R0,1);

Rmax_Matrix   = zeros(L_R0,L_tau);
Tosc_Matrix   = zeros(L_R0,L_tau);
R_R0_Matrix   = zeros(L_R0,L_tau);
Pbd_Matrix    = zeros(L_R0,L_tau);
Ubd_Matrix    = zeros(L_R0,L_tau);
tau_Matrix    = zeros(L_R0,L_tau);

%% Loop over array of R0
for i = 1:L_R0
    R0 = R0_Array(i);
    disp(['Simulating R0 = ' num2str(R0,'%.1e')])
    [tau_th,Rmax_th,Tosc_th, Results_Matrix] = inertial_conf_loop(R0,Rnbd_R0_ratio,tau_Array,Inertial_threshold);
    
    % Store the results
    tau_th_Array(i)  = tau_th;
    Rmax_th_Array(i) = Rmax_th;
    Tosc_th_Array(i) = Tosc_th;
    Rmax_Matrix(i,:) = Results_Matrix(1,:);
    Tosc_Matrix(i,:) = Results_Matrix(2,:);
    R_R0_Matrix(i,:) = Results_Matrix(3,:);
    Pbd_Matrix(i,:)  = Results_Matrix(4,:);
    Ubd_Matrix(i,:)  = Results_Matrix(5,:);
    tau_Matrix(i,:)  = Results_Matrix(6,:);
end


%% Plots, reproducing Figs. 3 and 4 in the Ref: Liang&Vogel2024 Preprint...
if L_R0 > 1 % contour map is only possible for array size of R0 >=2
fontsize = 14;
figure()
contourf(log10(tau_Matrix),log10(Rmax_Matrix),R_R0_Matrix)
xlabel('Log(\tau_L)')
ylabel('Log(R_{max})')
set(gca,'fontsize',fontsize)
title('Inertial confinement map')
c=colorbar;
c.Label.String = 'R|_{t=2\tau}/R_0';

figure()
loglog(tau_th_Array,Rmax_th_Array,'LineWidth',1.5)
xlabel('\tau_{L}')
ylabel('R_{max} (m)')
set(gca,'fontsize',fontsize)
title('Border of inert. conf. in (\tau_L,R_{max}) space')

figure()
loglog(Rmax_th_Array,Tosc_th_Array,'LineWidth',1.5)
xlabel('R_{max} (m)')
ylabel('T_{osc}')
set(gca,'fontsize',fontsize)
title('Border of inert. conf. in (R_{max},T_{osc}) space')

figure()
loglog(Rmax_th_Array,Tosc_th_Array./tau_th_Array,'LineWidth',1.5)
xlabel('R_{max} ')
ylabel('T_{osc}/\tau_L')
set(gca,'fontsize',fontsize)
end

%% Results for plotting in OriginLab; results by columns
tau_Matrix_Origin_col  = tau_Matrix(:);
Rmax_Matrix_Origin_col = Rmax_Matrix(:);
R_R0_Matrix_Origin_col = R_R0_Matrix(:);

% results by rows
tau_Matrix_Origin_row  = reshape(tau_Matrix.',1,[]);        % Originlab plot 2d contour map
Rmax_Matrix_Origin_row = reshape(Rmax_Matrix.',1,[]);
R_R0_Matrix_Origin_row = reshape(R_R0_Matrix.',1,[]);

log_tau_Matrix_Origin_row  = log10(tau_Matrix_Origin_row);  % Originlab plot 2d contour map
log_Rmax_Matrix_Origin_row = log10(Rmax_Matrix_Origin_row); % Originlab plot 2d contour map
