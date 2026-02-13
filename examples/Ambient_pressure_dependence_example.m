close all

p0 = logspace(5,7.69897,50); % ambient pressure 1e5 pa to 50e6 pa
P  = 1e9; % 1GPa pressure at bubble wall
fs = 14;  % fontsize

N  = length(p0);
Results = zeros(N,4); 

% Loop
for i=1:N
[rho0,c0,H,C] = TaitEOS(p0(i),'P',P);
Results(i,1)  = rho0;
Results(i,2)  = c0;
Results(i,3)  = H;
Results(i,4)  = C;
end

% Plots
figure()
semilogx(p0,Results(:,1),'LineWidth',1.5);
xlabel('Ambient Pressure (Pa)')
ylabel('Mass density of water \rho_{\infty}(kg/m^3)')
set(gca,'fontsize',fs)

figure()
semilogx(p0,Results(:,2),'LineWidth',1.5);
xlabel('Ambient Pressure (Pa)')
ylabel('Sound velocity c_{\infty} (m/s)')
set(gca,'fontsize',fs)

figure()
semilogx(p0,Results(:,3)/1e3,'LineWidth',1.5);
xlabel('Ambient Pressure (Pa)')
ylabel('Enthalpy \itH (kJ)')
legend(['Bubble wall pressure P = ' num2str(P/1e6) 'MPa'])
set(gca,'fontsize',fs)


figure()
semilogx(p0,Results(:,4),'LineWidth',1.5);
xlabel('Ambient Pressure (Pa)')
ylabel('Sound velocity at bubble wall \itC')
legend(['Bubble wall pressure \itP = ' num2str(P/1e6) 'MPa'])
set(gca,'fontsize',fs)
ylim([2739.5 2740.5])