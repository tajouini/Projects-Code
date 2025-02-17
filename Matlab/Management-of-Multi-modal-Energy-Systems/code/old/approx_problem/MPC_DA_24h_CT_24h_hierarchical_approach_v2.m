%%

%
% This code runs the predictive multi-energy management system.

%
% This upper level MPC represents the continuous trading planning finding place
% shortly before the delivery and closing at 12:00. The numerical values are taken form the following paper
%
%   "Multi-stage optimal energy amangement of multi-energy microgrid in
%    deregulated electricity markets", Wang et al., Applied
%               Energy(2022). 
%
% Copyright Institute of Automatic Control, Leibniz University Hannover.

%%
clear 
clc
close all
yalmip('clear')
addpath('c:\gurobi1002\win64\matlab')

%% Upper level MPC: Day ahead 
    Delta_DA = 1;                        % sampling time in [h]
    T_DA = 24;
    N_DA = T_DA/Delta_DA;                % discrete-time prediction horizon
    N24 = 24;
    N = N24;
 %% Lower level MPC: Continuous trading 
    Delta_CT = 1/4;                        % sampling time in [h]
    T_CT = 24; 
    N_CT = T_CT/Delta_CT;        % discrete-time prediction horizon [h]
    simStart = 1;
    sim_time = 24*4;               % in [h]
    simEnd = sim_time/Delta_CT;   % number of iterations   

    % size of the input vector u = [Pc; Pd; Hc; Hd]
    m = 4;
    % size of the state vector x = [SOC; H]
    n = 2;

 %% Setup figures
    f1 = figure(1);
    set(findall(gcf,'-property','FontSize'),'FontSize',20)
    set(gcf, 'Position',  [200, -90, 1000, 800])
    f2 = figure(2);
    set(findall(gcf,'-property','FontSize'),'FontSize',20)
    set(gcf, 'Position',  [200, -90, 1000, 800])
    sgtitle(f1,'Day ahead simulations');
    sgtitle(f2,'Closed-loop simulations');


   %% load measurement data
    % Day ahead (DA) measurement/forecast data
  % Load measurement data from 2010-01-01
    load ('C:\Users\Jouini\Documents\Research\IFES-Multi-modale-Energy-systeme\multi-energy-mangement-mpc\data\8760h-01p\PL_DA.mat');
    % Power measurements of PV generation 
    load('C:\Users\Jouini\Documents\Research\IFES-Multi-modale-Energy-systeme\multi-energy-mangement-mpc\data\8760h-01p\PV_DA.mat');
    % Power measurements of WT generation 
    load('C:\Users\Jouini\Documents\Research\IFES-Multi-modale-Energy-systeme\multi-energy-mangement-mpc\data\8760h-01p\WT_DA.mat');
    % Heat demand 
    load('C:\Users\Jouini\Documents\Research\IFES-Multi-modale-Energy-systeme\multi-energy-mangement-mpc\data\8760h-01p\HL_DA.mat');
    % Purchasing price (forecast)
    load('C:\Users\Jouini\Documents\Research\IFES-Multi-modale-Energy-systeme\multi-energy-mangement-mpc\data\8760h-01p\v_DA.mat'); 
    % Real purchasing price 
    load('C:\Users\Jouini\Documents\Research\IFES-Multi-modale-Energy-systeme\multi-energy-mangement-mpc\data\8760h-01p\v_DA_real.mat');
    % Selling price (see paper)
    u_DA = 0.9*v_DA;                                         
    % purchasing price of natural gas (see paper)
    v_GAS_DA = 0.055;  
    u_DA_real = 0.9*v_DA_real;
        



% Continuous Trading (CT) Read measurement/forecast data 
    % Electrical load measurement data 
    load ('C:\Users\Jouini\Documents\Research\IFES-Multi-modale-Energy-systeme\multi-energy-mangement-mpc\data\8760h-01p\PL_CT.mat');
    % Real load profile 
    %load ('C:\Users\Jouini\Documents\Research\IFES-Multi-modale-Energy-systeme\multi-energy-mangement-mpc\data\8760h-01p\PL_CT.mat');
    % Power measurements of PV generation 
    load('C:\Users\Jouini\Documents\Research\IFES-Multi-modale-Energy-systeme\multi-energy-mangement-mpc\data\8760h-01p\PV_CT.mat');
    % Power measurements of WT generation 
    load('C:\Users\Jouini\Documents\Research\IFES-Multi-modale-Energy-systeme\multi-energy-mangement-mpc\data\8760h-01p\WT_CT.mat');
    % Heat demand (TO CHECK) - I will just rescale Pload for now
    load('C:\Users\Jouini\Documents\Research\IFES-Multi-modale-Energy-systeme\multi-energy-mangement-mpc\data\8760h-01p\HL_CT.mat');
    % Purchasing price (forecast)
    load('C:\Users\Jouini\Documents\Research\IFES-Multi-modale-Energy-systeme\multi-energy-mangement-mpc\data\8760h-01p\v_CT.mat');
    % Real purchasing price
    load('C:\Users\Jouini\Documents\Research\IFES-Multi-modale-Energy-systeme\multi-energy-mangement-mpc\data\8760h-01p\v_CT_real.mat');
    % Selling price (see paper)
    u_CT = 0.9*v_CT;    
    % Purchasing price of natural gas (see paper)
    v_GAS_CT = 0.055; 
    % real puchasing price    
    u_CT_real = 0.8*v_CT_real;
    v_CT_real = 2*v_CT_real;



    %% Plot load and generation profiles
    % uncertainties_DA = [PL_DA(times); HL_DA(times); PV_DA(times); WT_DA(times); u_DA(times); v_DA(times)];
    figure(3)

    sgtitle('Uncertainty profiles for CT')
    set(findall(gcf,'-property','FontSize'),'FontSize',20)
    set(gcf, 'Position',  [100, -90, 1000, 800])
    
    ax31 = subplot(4,1,1);
    plot(v_DA(1:simEnd+N_CT))
    hold on     
    plot(u_DA(1:simEnd+N_CT))
    grid on
    title('Prices')
    xlabel('Step in (every 15 min)')
    ylabel('Price in dollars')
    xlim([1,simEnd+N_CT])

    ax32 = subplot(4,1,2);
    plot(PL_DA(1:simEnd+N_CT));
    grid on
    hold on
    plot(HL_DA(1:simEnd+N_CT));
    title('Load profile')
    xlabel('Step in (every 15 min)')
    ylabel('Power in kW')
    xlim([1,simEnd+N_CT])
    
    ax33 = subplot(4,1,3);
    plot(PV_DA(1:simEnd+N_CT));
    grid on
    title('PV generation')
    xlabel('Step in (every 15 min)')
    ylabel('Power in kW')
    xlim([1,simEnd+N_CT])

    ax34 = subplot(4,1,4);
    plot(WT_DA(1:simEnd+N_CT));
    grid on
    title('Wind Turbine generation')
    xlabel('Step in (every 15 min)')
    ylabel('Power in kW')
    xlim([1,simEnd+N_CT])




   
%% MPC Closed loop simulations
    [var_DA, state_DA] = generate_vars_v2(N_DA);
    [var_CT, state_CT] = generate_vars_v2(N_CT);
    % DA Initial conditions
    x_measure = [0.3;600];

    t = [];
    index_DA = 0;
    index_CT = 0;
    index_mon_DA = 0;
    index_mon_CT = 0;
    x = x_measure;
    u = [];
    Pgrid = [];
    Pg_DA_Cls = [];
    x_OL_vec = [];
    u_OL_vec = [];
    Pg_DA_vec = [];
    Pgrid_Cls = [];
    P_CHP_Cls = [];
    H_GB_Cls = [];
    rr = 1;
    s = 0;
   % Receding horizon control
    tic
for ii = simStart:simEnd
    clc
    %% Step 1: Day ahead
   if ii == rr % every 24h iterations, the upper level is triggered
    % update uncertainty profiles from the forecast for the next N=24h hours
    kk = 0;
    s = s+1;  % number of times DA is getriggered
    SoC_init_DA = x_measure(1);
    H_init_DA = x_measure(2);
    times = ceil(ii/4):N_DA-1+ceil(ii/4);  % take hourly measurements
    uncertainties_DA = [PL_DA(times); HL_DA(times); PV_DA(times); ...
        WT_DA(times); u_DA(times); v_DA(times)];
    % setup mpc for day ahead
    mpc_DA = setup_day_ahead_v2(var_DA, state_DA, uncertainties_DA, v_GAS_DA, N_DA, Delta_DA);
    % re-Initialisation 
    init_DA = {SoC_init_DA, H_init_DA};
    [sol_DA,flag_DA] = mpc_DA(init_DA);
    % out = {SoC,H,Pc,Pd,Hc,Hd,P_EB,P_CHP,I_ESS,I_HSS,H_EB,H_GB,Pgrid,cost_DA};
    if flag_DA > 0
        fprintf('Exitflag DA is %3d',flag_DA)
       break
    end
    % states
    
    x_OL(1,:) = sol_DA{1}; % SOC
    x_OL(2,:) = sol_DA{2}; % H
    % inputs
    u_OL(1,:) = sol_DA{3}; % Pc
    u_OL(2,:) = sol_DA{4}; % Pd
    u_OL(3,:) = sol_DA{5}; % Hc
    u_OL(4,:) = sol_DA{6}; % Hd
    Pg_DA = sol_DA{13}; % Pg

    P_DA_var = [sol_DA{3}; sol_DA{4}; Pg_DA; sol_DA{8}];
    H_DA_var = [sol_DA{5}; sol_DA{6}; sol_DA{12}];

    % DA performance index
    index_DA = index_DA + setup_cost_DA(v_GAS_DA,v_DA_real(times),u_DA_real(times),P_DA_var,H_DA_var, N24, Delta_DA);
    index_mon_DA = index_mon_DA + ...
    monetary_index_DA(v_GAS_DA,v_DA_real(times),u_DA_real(times),[sol_DA{13}; sol_DA{8}; sol_DA{12}], Delta_DA); 
    
    % interpolate DA solution from every hour to every 15 min
        % x_OL_intp = [];
        % for zz = 1:size(x_OL,2)-1
        %     x_OL_intp = [x_OL_intp interpolate_pts(x_OL(:,zz),u_OL(:,zz),Delta_DA,Delta_CT)];
        % end
        % x_OL_intp = [x_OL_intp x_OL(:,end)];

        % fixing the input trajectory to give it to the lower level
        %u_OL_rshp = [];
        Pg_DA_rshp = [];
        for jj = 1:size(u_OL,2)
          %u_OL_rshp = [u_OL_rshp repmat(u_OL(:,jj),1,Delta_DA*1/Delta_CT)]; 
          Pg_DA_rshp = [Pg_DA_rshp repmat(Pg_DA(jj),1,Delta_DA*1/Delta_CT)];
        end
        %u_OL_rshp = [u_OL_rshp u_OL(:,end) zeros(m,Delta_DA*24/Delta_CT)];
        Pg_DA_rshp = [Pg_DA_rshp Pg_DA(:,end)*ones(1,Delta_DA*1/Delta_CT)];
        rr = rr + 96;  % trigger DA every 24h --> 24*4 (in min)
        
        x_OL_vec = [x_OL_vec, x_OL];
        u_OL_vec = [u_OL_vec, u_OL];
        Pg_DA_vec = [Pg_DA_vec, Pg_DA];
      

       
        % plots of the Day ahead
        figure(1)

        subplot(2,4,1), cla
        hold on
        plot(1:N_DA*s+s, x_OL_vec(1,:),'b')
        grid on
        ylabel('SoC (in p.u.)')
        xlabel('Steps (in every 1h)')
        xlim([1, simEnd*Delta_CT])    
        drawnow

        subplot(2,4,2), cla
        hold on
        plot(1:N_DA*s+s, x_OL_vec(2,:),'r')
        grid on
        ylabel('H (in kWh)')
        xlabel('Steps (in every 1h)')
        xlim([1, simEnd*Delta_CT])    
        drawnow

        subplot(2,4,3), cla
        hold on
        plot(1:N_DA*s, u_OL_vec(1,:), 'b')
        grid on
        ylabel('Pc (in kW)')
        xlabel('Steps (in every 1h)')
        xlim([1, simEnd*Delta_CT])    
        drawnow

        subplot(2,4,4), cla
        hold on 
        plot(1:N_DA*s, u_OL_vec(2,:), 'b')
        grid on
        ylabel('Pd (in kW)')
        xlabel('Steps (in every 1h)')
        xlim([1, simEnd*Delta_CT])    
        drawnow

        subplot(2,4,5), cla
        hold on
        plot(1:N_DA*s, u_OL_vec(3,:),'r')
        grid on
        ylabel('Hc (in kW)')
        xlabel('Steps (in every 1h)')
        xlim([1, simEnd*Delta_CT])    
        drawnow

        subplot(2,4,6), cla
        hold on 
        plot(1:N_DA*s, u_OL_vec(4,:),'r')
        grid on
        ylabel('Hd (in kW)')
        xlabel('Steps (in every 1h)')
        xlim([1, simEnd*Delta_CT])    
        drawnow
 
        subplot(2,4,7), cla
        hold on
        plot(1:N_DA*s, Pg_DA_vec, 'g')
        hold on
        grid on
        ylabel('Pgrid (in kW)')
        xlabel('Steps (in every 1h)')
        xlim([1, simEnd*Delta_CT])    
        drawnow

    end
     %% step 2: Continuous trading
        SoC_init_CT = x_measure(1);
        H_init_CT = x_measure(2);
        % update uncertainty profiles from the CT forecast 
        times = ii:N_CT-1+ii;
        kk = kk+1;
        uncertainties_CT = [PL_CT(times); HL_CT(times); PV_CT(times); ...
        WT_CT(times); u_CT(times); v_CT(times)];
        % step 3: solve CT mpc problem
        mpc_CT = setup_continuous_trading_v2(var_CT,state_CT, uncertainties_CT, Pg_DA_rshp(kk:N_CT-1+kk), v_GAS_CT,N_CT, Delta_CT);
        % re-Initialisation 
        init_CT = {SoC_init_CT, H_init_CT};
        [sol_CT,flag_CT] = mpc_CT(init_CT);
        if flag_CT > 0
            fprintf('Exitflag CT is %3d',flag)
            break
        end
        
        % step 4: update the measurements
        % out = {SoC,H,Pc,Pd,Hc,Hd,P_EB,P_CHP,I_ESS,I_HSS,H_EB,H_GB(12),Pgrid(13),cost_CT};
        input = [sol_CT{3}(1);sol_CT{4}(1); sol_CT{5}(1); sol_CT{6}(1)];  % only first input in the sequence
        P_CT_var = [sol_CT{3}; sol_CT{4}; sol_CT{13}; sol_CT{8}];
        H_CT_var = [sol_CT{5}; sol_CT{6}; sol_CT{12}];
        [out1, out2] = dynamics(x_measure, input, Delta_CT); % simulate the closed-loop system
        x_measure = [out1; out2];
        x = [x, x_measure];
        u = [u, input];
        t = [t, Delta_CT];
        Pgrid_Cls = [Pgrid_Cls, sol_CT{13}(1)];
        P_CHP_Cls = [P_CHP_Cls, sol_CT{8}(1)];
        % H_GB_Cls = [H_GB_Cls, sol_CT{12}(1)];
        % Pg_DA_Cls = [Pg_DA_Cls, Pg_DA_rshp(1)];

         % CT performance index setup_cost_CT(v_GAS_CT,v,u,P_var,H_var,Pg_DA_rshp, N_CT, Delta_CT) 
        index_CT = index_CT + setup_cost_CT(v_GAS_CT,v_CT_real(ii),u_CT_real(ii),P_CT_var,H_CT_var, Pg_DA_rshp(kk), 1, Delta_CT);
        index_mon_CT = index_mon_CT + ...
        monetary_index_CT(v_GAS_CT,v_CT_real(ii),u_CT_real(ii),[sol_CT{13}(1); sol_CT{8}(1); sol_CT{12}(1)], Pg_DA_rshp(kk), Delta_CT); 

        % plots
        figure(2)
        subplot(2,4,1), cla
        hold on
        plot(x(1,:),'bo-')
        grid on
        ylabel('SoC (in p.u.)')
        xlabel('Steps (in every 15 min)')
        drawnow

        subplot(2,4,2), cla
        hold on
        plot(x(2,:),'ro-')
        grid on
        ylabel('H (in kWh)')
        xlabel('Steps (in every 15 min)')
        drawnow

        subplot(2,4,3), cla
        hold on
        plot(u(1,:),'bo-')
        grid on
        ylabel('Pc (in kW)')
        xlabel('Steps (in every 15 min)')
        drawnow

        subplot(2,4,4), cla
        hold on 
        plot(u(2,:),'bo-')
        grid on
        ylabel('Pd (in kW)')
        xlabel('Steps (in every 15 min)')
        drawnow

        subplot(2,4,5), cla
        hold on
        plot(u(3,:),'ro-')
        grid on
        ylabel('Hc (in kW)')
        xlabel('Steps (in every 15 min)')
        drawnow

        subplot(2,4,6), cla
        hold on 
        plot(u(4,:),'ro-')
        grid on
        ylabel('Hd (in kW)')
        xlabel('Steps (in every 15 min)')
        drawnow

        subplot(2,4,7), cla
        hold on
        plot(Pgrid_Cls,'go-')
        grid on
        ylabel('Pgrid (in kW)')
        xlabel('Steps (in every 15 min)')
        drawnow

        subplot(2,4,8), cla
        hold on
        plot(P_CHP_Cls,'go-')
        grid on
        ylabel('P_{CHP} (in kW)')
        xlabel('Steps (in every 15 min)')
        drawnow

        



 end
 toc



%% Print out the performance index
fprintf('Day ahead index is given by %3d and Continuous trading index is given by %3d \n', index_DA, index_CT);
fprintf('--------------------------------------------------------------------\n');
fprintf('Day ahead monetary index is given by %3d and continuous trading monetary index is given by %3d \n', index_mon_DA, index_mon_CT);
