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
    T_DA = 24;                           % operation time of the day ahead 
    T_SC = 12;                           % scheduling time of the day ahead
    N_pred_DA = T_DA/Delta_DA;             
    N_SC = T_SC/Delta_DA;                 
    N_DA = N_pred_DA + N_SC;	          % discrete-time prediction horizon 
    N24 = 24;
    N = N24;
 %% Lower level MPC: Continuous trading
    Delta_CT = 1/4;                        % sampling time in [h]
    T_CT = 1;
    N_CT = T_CT/Delta_CT;        % discrete-time prediction horizon [h]
    simStart = 1;
    days = 4;
    sim_time = 24*days;               % simulate for x days in [h]
    simEnd = sim_time/Delta_CT;   % number of iterations

    % size of the input vector u = [Pc; Pd; Hc; Hd]
    m = 4;
    % size of the state vector x = [SOC; H]
    n = 2;

 %% Setup figures
    f1 = figure(1);
    set(findall(gcf,'-property','FontSize'),'FontSize',20)
    sgtitle(f1,'Forecast uncertainty profiles for DA')
    set(findall(gcf,'-property','FontSize'),'FontSize',20)
    set(gcf, 'Position',  [100, -90, 500, 500])

    f2 = figure(2);
    set(findall(gcf,'-property','FontSize'),'FontSize',20)
    set(gcf, 'Position',  [200, -90, 1000, 800])
    sgtitle(f2,'24h-1h Day ahead simulations (next day)');
  

    %f3 = figure(3);
    %set(findall(gcf,'-property','FontSize'),'FontSize',20)
    %sgtitle(f3,'Final dispatch');
    %set(gcf, 'Position',  [100, -90, 1000, 800])

    f3 = figure(3);
    set(findall(gcf,'-property','FontSize'),'FontSize',20)
    sgtitle(f3,'24h-1h Closed-loop simulations (next day)');
    set(gcf, 'Position',  [100, -90, 1000, 800])


    noise = 0;
    p = load_parameters();


   %% load measurement data
load ('C:\Users\Jouini\Documents\Research\IFES-Multi-modale-Energy-systeme\multi-energy-mangement-mpc\data\8760h-01p\DA_uncertainties.mat');
load ('C:\Users\Jouini\Documents\Research\IFES-Multi-modale-Energy-systeme\multi-energy-mangement-mpc\data\8760h-01p\CT_uncertainties.mat');
load ('C:\Users\Jouini\Documents\Research\IFES-Multi-modale-Energy-systeme\multi-energy-mangement-mpc\data\8760h-01p\CT_real.mat');
load ('C:\Users\Jouini\Documents\Research\IFES-Multi-modale-Energy-systeme\multi-energy-mangement-mpc\data\8760h-01p\DA_real.mat');
 
    % v_DA = v_DA_N24;
    % PL_DA = PL_DA_N24;
    % HL_DA = HL_DA_N24;
    % PV_DA = PV_DA_N24;
    % WT_DA = WT_DA_N24;
    % 
    % v_CT = v_CT_N4;
    % PL_CT = PL_CT_N4;
    % HL_CT = HL_CT_N4;
    % PV_CT = PV_CT_N4;
    % WT_CT = WT_CT_N4;
  
    %% Plot load and generation profiles
    % uncertainties_DA = [PL_DA(times); HL_DA(times); PV_DA(times); WT_DA(times); u_DA(times); v_DA(times)];
    figure(1)
    
    ax31 = subplot(4,1,1);
    plot(v_DA(1:simEnd/4+N_DA))
    hold on
    plot(u_DA(1:simEnd/4+N_DA))
    grid on
    title('Prices')
    xlabel('Step in (every 1h)')
    ylabel('Price in dollars')
    xlim([1,simEnd/4+N_DA])

    ax32 = subplot(4,1,2);
    plot(PL_DA(1:simEnd/4+N_DA));
    grid on
    hold on
    plot(HL_DA(1:simEnd/4+N_DA));
    title('Load profile')
    xlabel('Step in (every 1h)')
    ylabel('Power in kW')
    xlim([1,simEnd/4+N_DA])

    ax33 = subplot(4,1,3);
    plot(PV_DA(1:simEnd/4+N_DA));
    grid on
    title('PV generation')
    xlabel('Step in (every 1h)')
    ylabel('Power in kW')
    xlim([1,simEnd/4+N_DA])

    ax34 = subplot(4,1,4);
    plot(WT_DA(1:simEnd/4+N_DA));
    grid on
    title('Wind Turbine generation')
    xlabel('Step in (every 1h)')
    ylabel('Power in kW')
    xlim([1,simEnd/4+N_DA])





%% MPC Closed loop simulations
    [var_DA, state_DA] = generate_vars(N_DA);
    [var_CT, state_CT] = generate_vars(N_CT);
    % DA Initial conditions measured 12h ago
    x_DA_measure = [0.3;600];
    x_CT_measure = [0.2592;100];

    t = [];
    index_DA = 0;
    index_CT = 0;
    index_mon_DA = 0;
    index_mon_CT = 0;
    x = [];
    u = [];
    Pgrid = [];
    Pg_DA_Cls = [];
    x_OL_vec = [];
    u_OL_vec = [];
    Pg_DA_vec = [];
    P_CHP_DA_vec = [];
    Pgrid_Cls = [];
    P_CHP_Cls = [];
    H_GB_Cls = [];
    
    P_CHP_final = [];
    Pc_final = [];
    Pd_final = [];
    Hc_final = [];
    Hd_final = [];
    rr = 1;
    zz = 1;
    s = 0;
   % Receding horizon control
    tic
for ii = simStart:simEnd
    clc
    %% Step I: Day ahead
    % assume we are at 12h one day before delivery, every 24h iterations, the upper level is triggered
   if ii == rr 
    % update uncertainty profiles from the forecast for the next N=24h hours
    kk = 0;
    s = s+1;  % number of times DA is getriggered
    SoC_init_DA = x_DA_measure(1);
    H_init_DA = x_DA_measure(2);
    times_DA = ceil(ii/4):N_DA-1+ceil(ii/4);  % take hourly measurements
    uncertainties_DA = [PL_DA(times_DA); HL_DA(times_DA); PV_DA(times_DA); ...
        WT_DA(times_DA); u_DA(times_DA); v_DA(times_DA)];
     %uncertainties_DA = [PL_DA(s,:); HL_DA(s,:); PV_DA(s,:); ...
     %  WT_DA(s,:); u_DA(s,:); v_DA(s,:)];
    % setup mpc for day ahead (var,state, uncertainties_DA, N_DA, p)
    mpc_DA = setup_day_ahead(var_DA, state_DA, uncertainties_DA, N_DA, p);
    % re-Initialisation
    init_DA = {SoC_init_DA, H_init_DA};
    [sol_DA,flag_DA] = mpc_DA(init_DA);
    % out = {SoC,H,Pc,Pd,Hc,Hd,P_EB,P_CHP,I_ESS,I_HSS,H_EB,H_GB,Pgrid,cost_DA};
    if flag_DA > 0
        fprintf('Exitflag DA is %3d',flag_DA)
       break
    end
    % states starting from the following day at midnight
    next_day = N_SC+1:N24+N_SC;
    x_OL(1,:) = sol_DA{1}(N_SC+1:end); % SOC
    x_OL(2,:) = sol_DA{2}(N_SC+1:end); % H
    % inputs
    u_OL(1,:) = sol_DA{3}(next_day); % Pc
    u_OL(2,:) = sol_DA{4}(next_day); % Pd
    u_OL(3,:) = sol_DA{5}(next_day); % Hc
    u_OL(4,:) = sol_DA{6}(next_day); % Hd
    Pg_DA = sol_DA{13}(next_day); % Pg
    P_CHP_DA = sol_DA{8}(next_day); % PCHP


    P_DA_var = [sol_DA{3}(next_day); sol_DA{4}(next_day); Pg_DA; sol_DA{8}(next_day)];
    H_DA_var = [sol_DA{5}(next_day); sol_DA{6}(next_day); sol_DA{12}(next_day)];
    mon_var = [sol_DA{13}(next_day); sol_DA{8}(next_day); sol_DA{12}(next_day)];
    % DA performance index
    times_P = ceil(ii/4)+N_SC:N_pred_DA-1+ceil(ii/4)+N_SC;  % start on the next day (i.e., after 12h)
    index_DA = index_DA + setup_cost_DA(v_DA_real(times_P),u_DA_real(times_P),P_DA_var,H_DA_var, N24, p);
    index_mon_DA = index_mon_DA + ...
    monetary_index_DA(v_DA_real(times_P),u_DA_real(times_P),mon_var, p);

    % interpolate DA solution from every hour to every 15 min (!)
        Pg_DA_rshp = [];
        for jj = 1:size(Pg_DA,2)
          Pg_DA_rshp = [Pg_DA_rshp repmat(Pg_DA(jj),1,Delta_DA*1/Delta_CT)];
        end
        Pg_DA_rshp = [Pg_DA_rshp Pg_DA_rshp(end)*ones(1,Delta_DA*1/Delta_CT)]; % for the horizon N_CT
        %Pg_DA_rshp = [Pg_DA_rshp  Pg_DA_rshp(1:N_CT)]; 
        
        rr = rr + 96;  % trigger DA every 24h --> 24*4 (in min)

        x_OL_vec = [x_OL_vec, x_OL];
        u_OL_vec = [u_OL_vec, u_OL];
        Pg_DA_vec = [Pg_DA_vec, Pg_DA];
        P_CHP_DA_vec = [P_CHP_DA_vec, P_CHP_DA];

        

      % plots of the Day ahead
        figure(2)

        subplot(2,4,1), cla
        hold on
        plot(1:N_pred_DA*s+s, x_OL_vec(1,:),'b')
        grid on
        ylabel('SoC (in p.u.)')
        xlabel('Steps (in every 1h)')
        xlim([1, simEnd*Delta_CT])
        drawnow

        subplot(2,4,2), cla
        hold on
        plot(1:N_pred_DA*s+s, x_OL_vec(2,:),'r')
        grid on
        ylabel('H (in kWh)')
        xlabel('Steps (in every 1h)')
        xlim([1, simEnd*Delta_CT])
        drawnow

        subplot(2,4,3), cla
        hold on
        plot(1:N_pred_DA*s, u_OL_vec(1,:), 'b')
        grid on
        ylabel('Pc (in kW)')
        xlabel('Steps (in every 1h)')
        xlim([1, simEnd*Delta_CT])
        drawnow

        subplot(2,4,4), cla
        hold on
        plot(1:N_pred_DA*s, u_OL_vec(2,:), 'b')
        grid on
        ylabel('Pd (in kW)')
        xlabel('Steps (in every 1h)')
        xlim([1, simEnd*Delta_CT])
        drawnow

        subplot(2,4,5), cla
        hold on
        plot(1:N_pred_DA*s, u_OL_vec(3,:),'r')
        grid on
        ylabel('Hc (in kW)')
        xlabel('Steps (in every 1h)')
        xlim([1, simEnd*Delta_CT])
        drawnow

        subplot(2,4,6), cla
        hold on
        plot(1:N_pred_DA*s, u_OL_vec(4,:),'r')
        grid on
        ylabel('Hd (in kW)')
        xlabel('Steps (in every 1h)')
        xlim([1, simEnd*Delta_CT])
        drawnow

        subplot(2,4,7), cla
        hold on
        plot(1:N_pred_DA*s, Pg_DA_vec, 'g')
        hold on
        grid on
        ylabel('Pgrid (in kW)')
        xlabel('Steps (in every 1h)')
        xlim([1, simEnd*Delta_CT])
        drawnow

        subplot(2,4,8), cla
        hold on
        plot(1:N_pred_DA*s, P_CHP_DA_vec, 'g')
        hold on
        grid on
        ylabel('P_{CHP} (in kW)')
        xlabel('Steps (in every 1h)')
        xlim([1, simEnd*Delta_CT])
        drawnow


    end
     %% Step II: Continuous trading
        % for the first time, take a closeby value of the predicted state
        % by DA
        % if ii==1   
        %    x_CT_measure = x_OL(1:2,1);
        % end
        SoC_init_CT = x_CT_measure(1);
        H_init_CT = x_CT_measure(2);
        % update uncertainty profiles from the CT forecast
        times = ii:N_CT-1+ii;
        kk = kk+1;
        %uncertainties_CT = [PL_CT(ii,:); HL_CT(ii,:); PV_CT(ii,:); ...
        %WT_CT(ii,:); u_CT(ii,:); v_CT(ii,:)];
        uncertainties_CT = [PL_CT(times); HL_CT(times); PV_CT(times); ...
        WT_CT(times); u_CT(times); v_CT(times)];
        % STEP 3: solve CT mpc problem
        mpc_CT = setup_continuous_trading(var_CT,state_CT, uncertainties_CT, Pg_DA_rshp(kk:N_CT-1+kk),N_CT,p);
        % re-Initialisation
        init_CT = {SoC_init_CT, H_init_CT};
        [sol_CT,flag_CT] = mpc_CT(init_CT);
        if flag_CT > 0
            fprintf('Exitflag CT is %3d',flag)
            break
        end

        % STEP 4: update the measurements
        % out = {SoC,H,Pc,Pd,Hc,Hd,P_EB,P_CHP,I_ESS,I_HSS,H_EB,H_GB(12),Pgrid(13),cost_CT};
        input = [sol_CT{3}(1);sol_CT{4}(1); sol_CT{5}(1); sol_CT{6}(1)];  % only first input in the sequence
        P_CT_var = [sol_CT{3}(1); sol_CT{4}(1); sol_CT{13}(1); sol_CT{8}(1)];
        H_CT_var = [sol_CT{5}(1); sol_CT{6}(1); sol_CT{12}(1)];
        vec = [sol_CT{13}(1); sol_CT{8}(1); sol_CT{12}(1)];
        %t = [t, Delta_CT];
        Pgrid_Cls = [Pgrid_Cls, sol_CT{13}(1)];
        P_CHP_Cls = [P_CHP_Cls, sol_CT{8}(1)];
        % H_GB_Cls = [H_GB_Cls, sol_CT{12}(1)];
        Pg_DA_Cls = [Pg_DA_Cls, Pg_DA_rshp(kk)];

        % CT performance index setup_cost_CT(v_GAS_CT,v,u,P_var,H_var,Pg_DA_rshp, N_CT, Delta_CT)
        index_CT = index_CT + setup_cost_CT(v_CT_real(ii),u_CT_real(ii),P_CT_var,H_CT_var, Pg_DA_rshp(kk), 1, p);
        index_mon_CT = index_mon_CT + ...
        monetary_index_CT(v_CT_real(ii),u_CT_real(ii),vec, Pg_DA_rshp(kk), p);

 
        [out1, out2] = dynamics(x_CT_measure, input, p); % simulate the closed-loop system
        x_CT_measure = [out1; out2];
        % update x_DA_measure every 12h
        zz = zz+ Delta_DA*12/Delta_CT;
        if ii==zz
            x_DA_measure = x_CT_measure;
        end
        x = [x, x_CT_measure];
        u = [u, input];
        


        % plots
        figure(3)
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
        hold on 
        plot(Pc_final,'ro-') 
        grid on
        ylabel('Pc (in kW)')
        xlabel('Steps (in every 15 min)')
        drawnow

        subplot(2,4,4), cla
        hold on
        plot(u(2,:),'bo-')
        hold on 
        plot(Pd_final,'ro-')
        grid on
        ylabel('Pd (in kW)')
        xlabel('Steps (in every 15 min)')
        drawnow

        subplot(2,4,5), cla
        hold on
        plot(u(3,:),'ro-')
        hold on
        plot(Hc_final,'bo-')
        grid on
        ylabel('Hc (in kW)')
        xlabel('Steps (in every 15 min)')
        drawnow

        subplot(2,4,6), cla
        hold on
        plot(u(4,:),'ro-')
        hold on 
        plot(Hd_final,'bo-')
        grid on
        ylabel('Hd (in kW)')
        xlabel('Steps (in every 15 min)')
        drawnow

        subplot(2,4,7), cla
        hold on
        plot(Pgrid_Cls,'bo-')
        hold on 
        plot(Pg_DA_Cls,'go-')
        grid on
        ylabel('Pgrid_{CT}/Pgrid_{DA} (in kW)' )
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
fprintf('Day ahead index is given by %3d and Continuous trading index is given by %3d\n', index_DA, index_CT);
fprintf('--------------------------------------------------------------------\n');
fprintf('Day ahead monetary index is given by %3d and continuous trading monetary index is given by %3d \n', index_mon_DA, index_mon_CT);

