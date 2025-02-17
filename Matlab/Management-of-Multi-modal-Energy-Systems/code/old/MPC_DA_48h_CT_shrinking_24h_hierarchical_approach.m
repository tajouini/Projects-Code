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
    T_DA = 48;
    N_DA = T_DA/Delta_DA;                % discrete-time prediction horizon
    N24 = 24;
    N = N24;
 %% Lower level MPC: Continuous trading
    Delta_CT = 1/4;                        % sampling time in [h]
    T_CT = 24;
    N_CT = T_CT/Delta_CT;        % discrete-time prediction horizon [h]
    simStart = 1;
    days = 4;
    sim_time = 24*days;
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
    sgtitle(f1,'48h-24h Day ahead simulations');
    sgtitle(f2,'48h-24h Closed-loop simulations');

    %% load measurement data
    % Day ahead (DA) measurement/forecast data
    % Load measurement 
    load ('C:\Users\Jouini\Documents\Research\IFES-Multi-modale-Energy-systeme\multi-energy-mangement-mpc\data\8760h-01p\PL_DA.mat');
    load ('C:\Users\Jouini\Documents\Research\IFES-Multi-modale-Energy-systeme\multi-energy-mangement-mpc\data\8760h-01p\PL_DA_N48.mat');
    % Power measurements of PV generation
    load('C:\Users\Jouini\Documents\Research\IFES-Multi-modale-Energy-systeme\multi-energy-mangement-mpc\data\8760h-01p\PV_DA.mat');
    load('C:\Users\Jouini\Documents\Research\IFES-Multi-modale-Energy-systeme\multi-energy-mangement-mpc\data\8760h-01p\PV_DA_N48.mat');
    % Power measurements of WT generation
    load('C:\Users\Jouini\Documents\Research\IFES-Multi-modale-Energy-systeme\multi-energy-mangement-mpc\data\8760h-01p\WT_DA.mat');
    load('C:\Users\Jouini\Documents\Research\IFES-Multi-modale-Energy-systeme\multi-energy-mangement-mpc\data\8760h-01p\WT_DA_N48.mat');
    % Heat demand
    load('C:\Users\Jouini\Documents\Research\IFES-Multi-modale-Energy-systeme\multi-energy-mangement-mpc\data\8760h-01p\HL_DA.mat');
    load('C:\Users\Jouini\Documents\Research\IFES-Multi-modale-Energy-systeme\multi-energy-mangement-mpc\data\8760h-01p\HL_DA_N48.mat');
    % Purchasing price (forecast)
    load('C:\Users\Jouini\Documents\Research\IFES-Multi-modale-Energy-systeme\multi-energy-mangement-mpc\data\8760h-01p\v_DA.mat');
    load('C:\Users\Jouini\Documents\Research\IFES-Multi-modale-Energy-systeme\multi-energy-mangement-mpc\data\8760h-01p\v_DA_N48.mat');
    % Real purchasing price
    load('C:\Users\Jouini\Documents\Research\IFES-Multi-modale-Energy-systeme\multi-energy-mangement-mpc\data\8760h-01p\v_DA_real.mat');
    
    
    


      % Continuous Trading (CT) Read measurement/forecast data
    % Electrical load measurement data
    load ('C:\Users\Jouini\Documents\Research\IFES-Multi-modale-Energy-systeme\multi-energy-mangement-mpc\data\8760h-01p\PL_CT.mat');
    load ('C:\Users\Jouini\Documents\Research\IFES-Multi-modale-Energy-systeme\multi-energy-mangement-mpc\data\8760h-01p\PL_CT_N96.mat');
    load('C:\Users\Jouini\Documents\Research\IFES-Multi-modale-Energy-systeme\multi-energy-mangement-mpc\data\8760h-01p\PL_CT_real.mat');
    % Power measurements of PV generation
    load('C:\Users\Jouini\Documents\Research\IFES-Multi-modale-Energy-systeme\multi-energy-mangement-mpc\data\8760h-01p\PV_CT.mat');
    load('C:\Users\Jouini\Documents\Research\IFES-Multi-modale-Energy-systeme\multi-energy-mangement-mpc\data\8760h-01p\PV_CT_N96.mat');
    load('C:\Users\Jouini\Documents\Research\IFES-Multi-modale-Energy-systeme\multi-energy-mangement-mpc\data\8760h-01p\PV_CT_real.mat');
    % Power measurements of WT generation
    load('C:\Users\Jouini\Documents\Research\IFES-Multi-modale-Energy-systeme\multi-energy-mangement-mpc\data\8760h-01p\WT_CT.mat');
    load('C:\Users\Jouini\Documents\Research\IFES-Multi-modale-Energy-systeme\multi-energy-mangement-mpc\data\8760h-01p\WT_CT_N96.mat');
    load('C:\Users\Jouini\Documents\Research\IFES-Multi-modale-Energy-systeme\multi-energy-mangement-mpc\data\8760h-01p\WT_CT_real.mat');
    % Heat demand (TO CHECK) - I will just rescale Pload for now
    load('C:\Users\Jouini\Documents\Research\IFES-Multi-modale-Energy-systeme\multi-energy-mangement-mpc\data\8760h-01p\HL_CT.mat');
    load('C:\Users\Jouini\Documents\Research\IFES-Multi-modale-Energy-systeme\multi-energy-mangement-mpc\data\8760h-01p\HL_CT_N96.mat');
    load('C:\Users\Jouini\Documents\Research\IFES-Multi-modale-Energy-systeme\multi-energy-mangement-mpc\data\8760h-01p\HL_CT_real.mat');
    % Purchasing price (forecast)
    load('C:\Users\Jouini\Documents\Research\IFES-Multi-modale-Energy-systeme\multi-energy-mangement-mpc\data\8760h-01p\v_CT.mat');
    load('C:\Users\Jouini\Documents\Research\IFES-Multi-modale-Energy-systeme\multi-energy-mangement-mpc\data\8760h-01p\v_CT_N96.mat');
    % Real purchasing price
    load('C:\Users\Jouini\Documents\Research\IFES-Multi-modale-Energy-systeme\multi-energy-mangement-mpc\data\8760h-01p\v_CT_real.mat');
    

    % v_DA = v_DA_N48;
    % PL_DA = PL_DA_N48;
    % HL_DA = HL_DA_N48;
    % PV_DA = PV_DA_N48;
    % WT_DA = WT_DA_N48;
    % 
    % v_CT = v_CT_N96;
    % PL_CT = PL_CT_N96;
    % HL_CT = HL_CT_N96;
    % PV_CT = PV_CT_N96;
    % WT_CT = WT_CT_N96;
    
    % Selling price (see paper)
    u_DA = 0.9*v_DA;
    % purchasing price of natural gas (see paper)
    v_GAS_DA = 0.055;
    u_DA_real = 0.9*v_DA_real;


     % Selling price (hier v_CT = interpolate(v_DA))
    u_CT = 0.8*0.9*v_CT;  % u_CT = 0.8*u_DA
    v_CT = 1.2*v_CT;  % v_CT = 1.2 v_DA
    % Purchasing price of natural gas (see paper)
    v_GAS_CT = 0.055;
    % Real purchasing price
    load('C:\Users\Jouini\Documents\Research\IFES-Multi-modale-Energy-systeme\multi-energy-mangement-mpc\data\8760h-01p\v_CT_real.mat');
    % real puchasing price
    u_CT_real = 0.8*0.9*v_CT_real;
    v_CT_real = 1.2*v_CT_real;


 %% load and generation profiles
    % uncertainties_DA = [PL_DA(times); HL_DA(times); PV_DA(times); WT_DA(times); u_DA(times); v_DA(times)];
    figure(3)

    sgtitle('Forecast uncertainty profile for DA')
    set(findall(gcf,'-property','FontSize'),'FontSize',20)
    set(gcf, 'Position',  [100, -90, 1000, 800])

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
    xlim([1,simEnd/4])

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

    % DA Initial conditions
    x_measure = [0.3;600];
    t = [];
    x = x_measure;
    u = [];
    Pgrid_Cls = [];
    P_CHP_Cls = [];
    H_GB_Cls = [];
    Pg_DA_Cls = [];

    x_OL_vec = [];
    u_OL_vec = [];
    Pg_DA_vec = [];
    rr = 1;
    s = 0;
    index_DA = 0;
    index_CT = 0;
    index_mon_DA = 0;
    index_mon_CT = 0;
   % Receding horizon control
    tic
for ii = simStart:simEnd   % in every 1/4 hour
    clc
    %disp(t)
    %% Step 1: Day ahead
   if ii == rr % every 24h iterations, the upper level is triggered
    kk = 0;
    s = s+1;  % number of times DA is getriggered
    SoC_init_DA = x_measure(1);
    H_init_DA = x_measure(2);
    times_DA = ceil(ii/4):N_DA-1+ceil(ii/4);  % take hourly measurements
    times24 = ceil(ii/4):N24-1+ceil(ii/4);
    % update uncertainty profiles from the forecast for the next N=24h hours
    % uncertainties_DA = [PL_DA(s,:); HL_DA(s,:); PV_DA(s,:); ...
    %     WT_DA(s,:); u_DA(s,:); v_DA(s,:)];
    uncertainties_DA = [PL_DA(times_DA); HL_DA(times_DA); PV_DA(times_DA); ...
        WT_DA(times_DA); u_DA(times_DA); v_DA(times_DA)];
    % setup mpc for day ahead
    mpc_DA = setup_day_ahead(var_DA, state_DA, uncertainties_DA, v_GAS_DA, N_DA, Delta_DA);
    % re-Initialisation
    init_DA = {SoC_init_DA, H_init_DA};
    [sol_DA,flag_DA] = mpc_DA(init_DA);
    if flag_DA > 0
        fprintf('Exitflag DA is %3d',flag_DA)
       break
    end
    % states
    x_OL(1,:) = sol_DA{1}(1:25); % take only first 24h predictions
    x_OL(2,:) = sol_DA{2}(1:25);
    %inputs
    u_OL(1,:) = sol_DA{3}(1:24); % take only first 24h predictions
    u_OL(2,:) = sol_DA{4}(1:24);
    u_OL(3,:) = sol_DA{5}(1:24);
    u_OL(4,:) = sol_DA{6}(1:24);
    Pg_DA = sol_DA{13}; % P_g


    P_DA_var = [sol_DA{3}(1:24); sol_DA{4}(1:24); sol_DA{13}(1:24); sol_DA{8}(1:24)];
    H_DA_var = [sol_DA{5}(1:24); sol_DA{6}(1:24); sol_DA{12}(1:24)];

    % DA performance index
    index_DA = index_DA + setup_cost_DA(v_GAS_DA,v_DA_real(times24),u_DA_real(times24),P_DA_var,H_DA_var, N24, Delta_DA);
    index_mon_DA = index_mon_DA + ...
    monetary_index_DA(v_GAS_DA,v_DA_real(times24),u_DA_real(times24),[sol_DA{13}; sol_DA{8}; sol_DA{12}], Delta_DA);

    % fixing the input trajectory to give it to the lower level
    Pg_DA_rshp = [];
    for jj = 1:size(Pg_DA,2)
       Pg_DA_rshp = [Pg_DA_rshp repmat(Pg_DA(jj),1,Delta_DA*1/Delta_CT)];
    end
       %u_OL_rshp = [u_OL_rshp u_OL(:,end)*ones(m,Delta_DA*24/Delta_CT)];
       %Pg_DA_rshp = [Pg_DA_rshp Pg_DA(:,end)*ones(1,Delta_DA*24/Delta_CT)];
       rr = rr + Delta_DA/Delta_CT*24;  % trigger DA every 24h --> 24*4 (in min)

       N_CT = T_CT/Delta_CT; % reset the horizon
       [var_CT, state_CT] = generate_vars(N_CT); % regenerate CT vars
        x_OL_vec = [x_OL_vec, x_OL];
        u_OL_vec = [u_OL_vec, u_OL];
        Pg_DA_vec = [Pg_DA_vec, Pg_DA];

        % plots of the Day ahead
        figure(1)

        subplot(2,4,1), cla
        hold on
        plot(1:N24*s+s, x_OL_vec(1,:),'b')
        grid on
        ylabel('SoC (in p.u.)')
        xlabel('Steps (in every 1h)')
        xlim([1, simEnd*Delta_CT])
        drawnow

        subplot(2,4,2), cla
        hold on
        plot(1:N24*s+s, x_OL_vec(2,:),'r')
        grid on
        ylabel('H (in kW)')
        xlabel('Steps (in every 1h)')
        xlim([1, simEnd*Delta_CT])
        drawnow

        subplot(2,4,3), cla
        hold on
        plot(1:N24*s, u_OL_vec(1,:), 'b')
        grid on
        ylabel('Pc (in kW)')
        xlabel('Steps (in every 1h)')
        xlim([1, simEnd*Delta_CT])
        drawnow

        subplot(2,4,4), cla
        hold on
        plot(1:N24*s, u_OL_vec(2,:), 'b')
        grid on
        ylabel('Pd (in kW)')
        xlabel('Steps (in every 1h)')
        xlim([1, simEnd*Delta_CT])
        drawnow

        subplot(2,4,5), cla
        hold on
        plot(1:N24*s, u_OL_vec(3,:),'r')
        grid on
        ylabel('Hc (in kW)')
        xlabel('Steps (in every 1h)')
        xlim([1, simEnd*Delta_CT])
        drawnow

        subplot(2,4,6), cla
        hold on
        plot(1:N24*s, u_OL_vec(4,:),'r')
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
        kk = kk+1;
        times = ii:N_CT-1+ii;
         if (kk+N_CT-1>length(Pg_DA_rshp)) %&& (N_CT>1)
         N_CT = N_CT-1;
         [var_CT, state_CT] = generate_vars(N_CT);
         end
        times = ii:N_CT-1+ii;
        % uncertainties_CT = [PL_CT(ii,1:N_CT); HL_CT(ii,1:N_CT); PV_CT(ii,1:N_CT); ...
        % WT_CT(ii,1:N_CT); u_CT(ii,1:N_CT); v_CT(ii,1:N_CT)];
        uncertainties_CT = [PL_CT(times); HL_CT(times); PV_CT(times); ...
        WT_CT(times); u_CT(times); v_CT(times)];
        % step 3: solve CT mpc problem
        mpc_CT = setup_continuous_trading(var_CT,state_CT, uncertainties_CT, Pg_DA_rshp(kk:kk+N_CT-1), v_GAS_CT,N_CT, Delta_CT);
        % re-Initialisation
        init_CT = {SoC_init_CT, H_init_CT};
        [sol_CT,flag_CT] = mpc_CT(init_CT);
        if flag_CT > 0
            fprintf('Exitflag CT is %3d',flag)
            break
        end

        % step 4: update the measurements
        % out = {SoC,H,Pc,Pd,Hc,Hd,P_EB,P_CHP,I_ESS,I_HSS,H_EB,H_GB,Pgrid,cost_CT};
        P_CT_var = [sol_CT{3}(1); sol_CT{4}(1); sol_CT{13}(1); sol_CT{8}(1)];
        H_CT_var = [sol_CT{5}(1); sol_CT{6}(1); sol_CT{12}(1)];
        input = [sol_CT{3}(1);sol_CT{4}(1); sol_CT{5}(1); sol_CT{6}(1)];  % only first input in the sequence
        [out1,out2] = dynamics(x_measure, input, Delta_CT);
        x_measure = [out1; out2];
        x = [x, x_measure];
        u = [u, input];
        t = [t, Delta_CT];
        Pgrid_Cls = [Pgrid_Cls, sol_CT{13}(1)];
        P_CHP_Cls = [P_CHP_Cls, sol_CT{8}(1)];
        %H_GB_Cls = [H_GB_Cls, sol_CT{12}(1)];
        Pg_DA_Cls = [Pg_DA_Cls, Pg_DA_rshp(kk)];

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
        ylabel('H (in kW)')
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
        hold on 
        plot(Pg_DA_Cls,'bo-')
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
% print out
fprintf('Day ahead index is given by %3d and continuous trading index is given by %3d', index_DA, index_CT);
fprintf('\n--------------------------------------------------------------------\n');
fprintf('Day ahead monetary index is given by %3d and continuous trading monetary index is given by %3d \n', index_mon_DA, index_mon_CT);



