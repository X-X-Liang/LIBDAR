close all
clearvars

%% Show test header (useful when run as part of test suite)
disp('-----------------------------------------------------------------')
disp('-----            LIBDAR TOOLBOX TEST:                       -----')
disp('-----      Bubble dynamics and acoustic radiation example   -----')
disp('-----   Ref: Tian et al. submitted to PRL... (2026) -----')
disp(' ')

%% Add path to necessary functions
toolboxPath=fileparts(fileparts(mfilename('fullpath'))); %Get the toolbox path

%Add core, misc and expData path
addpath(fullfile(toolboxPath,'core_functions'),fullfile(toolboxPath,'misc'),fullfile(toolboxPath,'expData'));

%% explicit input parameters
tp   = 6e-9;      % effective pulse duration
p0   = 50e6;       % ambient pressure

disp('Simulating...')
tic

%% load experimental data
expdata = xlsread('expData_0.1MPa-50MPa.xlsx');
switch p0
    case 1e5
        time_exp   = expdata(:,1)/1e3; % in us
        radius_exp = expdata(:,2);     % in um
        T1   = 201.147e-6; % measured 1st osc. time
        T2   = 63.097e-6; % measured 2nd osc. time
        R0   = 20e-6;     % initial bubble radius, taken as exp. determined plasma radius
        Rnbd = 400.4e-6;    % equilibrium bubble radius after bd, tuning parameter such that simulated Tosc1 fits with measured value
        Rnc1 = 132e-6;    % equilibrium bubble radius at 1st bubble collapse, tuning parameter such that simulated Tosc2 fits with measured value
        Rnc2 = 132e-6;    % equilibrium bubble radius at 2nd bubble collapse, tuning parameter such that simulated Tosc3 fits with measured value

    case 2e6
        time_exp   = expdata(:,3)/1e3; % in us
        radius_exp = expdata(:,4);     % in um
        T1   = 14.409e-6; 
        T2   = 6.078e-6;  
        R0   = 18e-6;     
        Rnbd = 166.5e-6;   
        Rnc1 = 69e-6;
        Rnc2 = 69e-6;

    case 5e6
        time_exp   = expdata(:,5)/1e3; % in us
        radius_exp = expdata(:,6);     % in um
        T1   = 6.901e-6;
        T2   = 3.189e-6;  
        R0   = 17e-6;     
        Rnbd = 136.1e-6;   
        Rnc1 = 61e-6;
        Rnc2 = 61e-6;

    case 10e6
        time_exp   = expdata(:,7)/1e3; % in us
        radius_exp = expdata(:,8);     % in um
        T1   = 4.002e-6;
        T2   = 2.273e-6;  
        R0   = 17e-6;     
        Rnbd = 117.7e-6;   
        Rnc1 = 65e-6;
        Rnc2 = 65e-6;

    case 15e6
        time_exp   = expdata(:,9)/1e3; % in us
        radius_exp = expdata(:,10);    % in um
        T1   = 2.902e-6;
        T2   = 1.857e-6;
        R0   = 17e-6;
        Rnbd = 107.3e-6;
        Rnc1 = 67e-6;
        Rnc2 = 67e-6;

    case 20e6
        time_exp   = expdata(:,11)/1e3; % in us
        radius_exp = expdata(:,12);     % in um
        T1   = 2.286e-6;
        T2   = 1.534e-6;
        R0   = 17e-6;
        Rnbd = 99.3e-6;
        Rnc1 = 65e-6;
        Rnc2 = 65e-6;

    case 25e6
        time_exp   = expdata(:,13)/1e3; % in us
        radius_exp = expdata(:,14);     % in um
        T1   = 1.9e-6;
        T2   = 1.393e-6;
        R0   = 17e-6;
        Rnbd = 92.8e-6;
        Rnc1 = 67e-6;
        Rnc2 = 67e-6;


    case 30e6
        time_exp   = expdata(:,15)/1e3; % in us
        radius_exp = expdata(:,16);     % in um
        T1   = 1.639e-6;
        T2   = 1.313e-6;
        R0   = 17e-6;
        Rnbd = 87.4e-6;
        Rnc1 = 70e-6;
        Rnc2 = 70e-6;

    case 35e6
        time_exp   = expdata(:,17)/1e3; % in us
        radius_exp = expdata(:,18);     % in um
        T1   = 1.456e-6;
        T2   = 1.208e-6;
        R0   = 17e-6;
        Rnbd = 83.9e-6;
        Rnc1 = 70e-6;
        Rnc2 = 70e-6;

    case 40e6
        time_exp   = expdata(:,19)/1e3; % in us
        radius_exp = expdata(:,20);     % in um
        T1   = 1.302e-6;
        T2   = 1.125e-6;
        R0   = 16.5e-6;
        Rnbd = 80e-6;
        Rnc1 = 70e-6;
        Rnc2 = 70e-6;

    case 45e6
        time_exp   = expdata(:,21)/1e3; % in us
        radius_exp = expdata(:,22);     % in um
        T1   = 1.192e-6;
        T2   = 1.059e-6;
        R0   = 16.5e-6;
        Rnbd = 77.5e-6;
        Rnc1 = 70e-6;
        Rnc2 = 70e-6;

    case 50e6
        time_exp   = expdata(:,23)/1e3; % in us
        radius_exp = expdata(:,24);     % in um
        T1   = 1.106e-6;
        T2   = 1.01e-6;
        R0   = 16.2e-6;
        Rnbd = 75.5e-6;
        Rnc1 = 70.5e-6;
        Rnc2 = 70.5e-6;

    otherwise
        msg = ['Check your input please! Experimental data not available at p_0 = ' num2str(p0/1e6) 'MPa'];
        error(msg);

end
%% bubble dynamics

[filename,t,R,P,U,tind_sw,tindRmin] = bubble_dynamics_and_energy_partitioning (T1,T2,R0,Rnbd,Rnc1,Rnc2,tp,'p0',p0,'GS','adaptive','EP',0,'JSC','general');

% exp. v.s. cal. curves
fontsize = 14;
figure();
plot(time_exp, radius_exp, 'ro',t*1e6,R*1e6,'b-','LineWidth',1.5);
xlabel('Time (\mus)')
ylabel('Radius (\mum)')
legend('exp.','theo')
title(['R(t) exp. v.s. theo. with p_0 = ' num2str(p0/1e6) ' MPa']);
set(gca,'fontsize',fontsize)

%% acoustic radiation

switch p0
    
    case {1e5, 5e6, 50e6}

    isARexp = 1; % is acoustic radiation at bubble expansion? 1 for bubble expansion and 0 for collapse and rebound

    [tindstart,tindend,dtkb,r_container,u_container,p_container,Rtt,Utt,Ptt] = KB_acoustic_radiation(t,R,P,U,tind_sw,'p0',p0,'isARexp',isARexp,'winSel','manuel','dtkb',2e-12,'t_sw0',t(2),'t_swd',200e-9);
 
    % acoustic radiation, animation, and publication quality plots, reproducing Figs. 9 and 10 in Liang2022JFM
    time_q_pub = [3:1:50,52:2:104];% querried time in unit ns

    [~,~,~,~,~,SW_radius_time_curve,SW_velocity_time_curve] = acoustRadwReshaping(t,R,tp,tindstart,tindend,tindRmin,dtkb,r_container,u_container,p_container,Rtt,Utt,Ptt,'time_q',time_q_pub,'isARexp',isARexp,'xylim',[10 1000 0.1 1e4],'Nsel',8,'isExtrafig',0);
    
    switch p0
        case 1e5
            sw_time_exp   = expdata(:,25); % in ns
            sw_radius_exp = expdata(:,26); % in um
        case 5e6
            sw_time_exp   = expdata(:,27); % in ns
            sw_radius_exp = expdata(:,28); % in um

        case 50e6
            sw_time_exp   = expdata(:,29); % in ns
            sw_radius_exp = expdata(:,30); % in um

    end

    figure()
    plot(sw_time_exp, sw_radius_exp, 'bo',SW_radius_time_curve(:,1),SW_radius_time_curve(:,2),'b','LineWidth',1.5)
    xlabel('Time (ns)','fontsize',fontsize)
    ylabel('Shock Wave Radius (μm)','fontsize',fontsize)
    legend('exp.','theo')
    set(gca,'fontsize',fontsize)

    figure()
    plot(SW_velocity_time_curve(:,1),SW_velocity_time_curve(:,2),'b','LineWidth',1.5)
    xlabel('Time (ns)','fontsize',fontsize)
    ylabel('Shock Wave Velocity (μm)','fontsize',fontsize)
    set(gca,'fontsize',fontsize)

end

toc
%% Add blank line (nicer formatting for test text output).
disp('-----                  End of simulation                -------- ')