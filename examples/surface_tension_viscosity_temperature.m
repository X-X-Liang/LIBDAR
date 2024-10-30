close all;
clearvars

%% Show test header (useful when run as part of test suite)
disp('-----------------------------------------------------------------')
disp('-----            LIBDAR TOOLBOX TEST:                       -----')
disp('-----          sigma_T and mu_T test                        -----')
disp('-----------------------------------------------------------------')
disp(' ')

%% Add path to necessary functions
toolboxPath=fileparts(fileparts(mfilename('fullpath'))); %Get the toolbox path

%Add core and misc path
addpath(fullfile(toolboxPath,'core_functions'),fullfile(toolboxPath,'misc'));

T = 273:647;   % array for temperature from 0 degree Celsius to the critical point
sigma = 0*T;   % empty array
mu    = 0*T;   % empty array
pv    = 0*T;   % empty array
fonts = 14;    % fontsize

for i=1:length(T)
    Temperature = T(i);
    sigma(i)    = sigma_T(Temperature); % surface tension
    mu(i)       = mu_T(Temperature);    % viscosity
    pv(i)       = pv_T(Temperature);    % saturated vapor pressre
end

%% plots
figure()
plot(T,sigma,'LineWidth',1.5)
xlabel('Temperature (K)');
ylabel('Surface tension (N/m)');
set(gca,'fontsize',fonts);

figure()
plot(T,mu,'LineWidth',1.5)
xlabel('Temperature (K)');
ylabel('Viscosity (Pa.s)');
set(gca,'fontsize',fonts);

figure()
semilogy(T,pv,'LineWidth',1.5)
xlabel('Temperature (K)');
ylabel('Pressure (Pa)');
set(gca,'fontsize',fonts);