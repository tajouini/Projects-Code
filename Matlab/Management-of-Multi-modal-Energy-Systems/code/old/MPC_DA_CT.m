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
function [index_DA, index_CT, index_mon_DA, index_mon_CT] = MPC_DA_CT (varargin)

close all
yalmip('clear')

p = struct ();

%% Input Arguments
    p.T_DA = 48;
    p.T_CT = 1;
    p.days = 4;
    if nargin > 0
      p.T_DA = varargin{1};
      if nargin > 1
        p.T_CT = varargin{2};
        if nargin > 2
          p.days = varargin{3};
        end
      end
    end

%% Upper level MPC: Day ahead
    p.Delta_DA = 1;                        % sampling time in [h]
    p.N_DA = p.T_DA/p.Delta_DA;                % discrete-time prediction horizon
    p.N24 = 24;
%% Lower level MPC: Continuous trading
    p.Delta_CT = 1/4;                        % sampling time in [h]
    p.N_CT = p.T_CT/p.Delta_CT;        % discrete-time prediction horizon [h]
    simStart = 1;
    sim_time = 24*p.days;               % simulate for a month in [h]
    simEnd = sim_time/p.Delta_CT;   % number of iterations

    % size of the input vector u = [Pc; Pd; Hc; Hd]
    m = 4;
    % size of the state vector x = [SOC; H]
    n = 2;

%% Other system constants
    p.Eta_EB = 0.95;    % Electricity to heat efficiency of the EB
    p.b = 1.2;          % Characteristic of the heat-to-electric ratio
    p.Pmax = 500;       % Maximum charging/discharging power
    p.Q_ESS = 2000;     % Capacity of ESS in [kWh]
    p.EtaC = 0.95;      % Charging efficiency in p.u.
    p.EtaD = 0.95;      % Discharging efficiency in p.u.
    p.Rho_ESS = 0.01;   % Depreciation coefficient in $/kWh
    p.Rho_HSS = 0.01;   % Depreciation coefficient in $/kWh
    p.Eta_CHP = 0.4;    %0.95; % Gas to electricity efficiency of CHP
    p.Eta_GB = 0.8;     % Gas to heat efficiency of the GB
    p.Eta_HC = 0.9;     % Heat charging efficiency in p.u.
    p.Eta_HD = 0.9;     % Discharging efficiency in p.u.
    p.SocMax = 0.8;     % Maximal SOC
    p.SocMin = 0.2;     % Minimal SOC
    p.Pmax_CHP = 1200;  % Maximal combined power to heat (in kW)
    p.Pmin_CHP = 0;     % minimal combined power to heat (in kW)
    p.HSS_Max = 3000;   % Maximal stored heat (in kWh)
    p.HSS_Min = 100;    % Minimal stored heat (in kWh)
    p.Hmax = 750;       % Maximum chargin/discharding power
    p.H_GBMax = 500;    % Maximal heat power limit of the GB
    p.H_GBMin = 0;      % Minimal heat power limit of the GB
    p.H_EBMax = 500;    % Maximal heat power limit of the EB
    p.H_EBMin = 0;      % Minimal heat power limit of the EB
    p.PgridMax = 2500;  % Maximal grid power
    p.PgridMin = -2500; % Minimal grid power
    p.v_GAS_DA = 0.055; % Purchasing price of natural gas (see paper)
    p.v_GAS_CT = 0.055;

%% Setup figures
    f1 = figure(1);
    set(findall(gcf,'-property','FontSize'),'FontSize',20)
    set(gcf, 'Position',  [200, -90, 1000, 800])
    f2 = figure(2);
    set(findall(gcf,'-property','FontSize'),'FontSize',20)
    set(gcf, 'Position',  [200, -90, 1000, 800])


   %% load measurement data
    % Day ahead (DA) measurement/forecast data
  % Load measurement data from 2010-01-01
    load ('../data/8760h-01p/PL_DA.mat');
    % Power measurements of PV generation
    load('../data/8760h-01p/PV_DA.mat');
    % Power measurements of WT generation
    load('../data/8760h-01p/WT_DA.mat');
    % Heat demand
    load('../data/8760h-01p/HL_DA.mat');
    % Purchasing price (forecast)
    load('../data/8760h-01p/v_DA.mat');
     % Real purchasing price
    load('../data/8760h-01p/v_DA_real.mat');




% Continuous Trading (CT) Read measurement/forecast data
    % Electrical load measurement data
    load ('../data/8760h-01p/PL_CT.mat');
    % Power measurements of PV generation
    load('../data/8760h-01p/PV_CT.mat');
    % Power measurements of WT generation
    load('../data/8760h-01p/WT_CT.mat');
    % Heat demand (TO CHECK) - I will just rescale Pload for now
    load('../data/8760h-01p/HL_CT.mat');
    % Purchasing price (forecast)
    load('../data/8760h-01p/v_CT.mat');
     % Real purchasing price
    load('../data/8760h-01p/v_CT_real.mat');

    % Selling price (see paper)
    u_DA = 0.9*v_DA;
    u_DA_real = 0.9*v_DA_real;

    % Selling price
    u_CT = 0.8*v_CT;  % see paper
    v_CT = 2*v_CT;  % see paper
    % real puchasing price (same as in the forecast)
    u_CT_real = 0.8*v_CT_real;
    v_CT_real = 2*v_CT_real;

    %% Plot load and generation profiles
    % uncertainties_DA = [PL_DA(times); HL_DA(times); PV_DA(times); WT_DA(times); u_DA(times); v_DA(times)];
    figure(3)
    set(findall(gcf,'-property','FontSize'),'FontSize',20)
    set(gcf, 'Position',  [100, -90, 1000, 800])

    ax31 = subplot(4,1,1);
    plot(v_DA(1:simEnd/4+p.N_DA))
    hold on
    plot(u_DA(1:simEnd/4+p.N_DA))
    grid on
    set( ax31, 'title', 'Prices' );
    xlabel('Step in (every 1h)')
    ylabel('Price in dollars')
    xlim([1,simEnd/4+p.N_DA])

    ax32 = subplot(4,1,2);
    plot(PL_DA(1:simEnd/4+p.N_DA));
    grid on
    hold on
    plot(HL_DA(1:simEnd/4+p.N_DA));
    set( ax32, 'title', 'Load profiles' );
    xlabel('Step in (every 1h)')
    ylabel('Power in kW')
    xlim([1,simEnd/4+p.N_DA])

    ax33 = subplot(4,1,3);
    plot(PV_DA(1:simEnd/4+p.N_DA));
    grid on
    set( ax33, 'title', 'PV generation' );
    xlabel('Step in (every 1h)')
    ylabel('Power in kW')
    xlim([1,simEnd/4+p.N_DA])

    ax34 = subplot(4,1,4);
    plot(WT_DA(1:simEnd/4+p.N_DA));
    grid on
    set( ax34, 'title', 'Wind turbine generation' );
    xlabel('Step in (every 1h)')
    ylabel('Power in kW')
    xlim([1,simEnd/4+p.N_DA])

    S3 = axes( 'visible', 'off', 'title', 'Forecast uncertainty profiles for DA', 'FontSize', 16);




%% MPC Closed loop simulations
    [var_DA, state_DA] = generate_vars(p.N_DA);
    [var_CT, state_CT] = generate_vars(p.N_CT);
    % DA Initial conditions
    x_measure = [0.3;600];

    t = [];

    index_DA = zeros (1,p.days);
    index_mon_DA = zeros (1,p.days);
    index_CT = 0;
    index_mon_CT = 0;

    x = x_measure;

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

    rr = 1;
    s = 0;
   % Receding horizon control
    tic

for ii = simStart:simEnd
    %% Step 1: Day ahead
   if ii == rr % every 24h iterations, the upper level is triggered
    % update uncertainty profiles from the forecast for the next N=24h hours
    kk = 0;
    s = s+1;  % number of times DA is getriggered
    SoC_init_DA = x_measure(1);
    H_init_DA = x_measure(2);
    times = ceil(ii/4):p.N_DA-1+ceil(ii/4);  % take hourly measurements
    times24 = ceil(ii/4):p.N24-1+ceil(ii/4);
    uncertainties_DA = [PL_DA(times); HL_DA(times); PV_DA(times); WT_DA(times); u_DA(times); v_DA(times)];
    % setup mpc for day ahead
    mpc_DA = setup_day_ahead(var_DA, state_DA, uncertainties_DA, p);
    % re-Initialisation
    init_DA = {SoC_init_DA, H_init_DA};
    [sol_DA,flag_DA] = mpc_DA(init_DA);
    % out = {SoC,H,Pc,Pd,Hc,Hd,P_EB,P_CHP,I_ESS,I_HSS,H_EB,H_GB,Pgrid,cost_DA};
    if flag_DA > 0
        fprintf('Exitflag DA is %3d',flag_DA)
       break
    end
    % states

    x_OL(1,:) = sol_DA{1}(1:25); % SOC
    x_OL(2,:) = sol_DA{2}(1:25); % H
    % inputs
    u_OL(1,:) = sol_DA{3}(1:24); % Pc
    u_OL(2,:) = sol_DA{4}(1:24); % Pd
    u_OL(3,:) = sol_DA{5}(1:24); % Hc
    u_OL(4,:) = sol_DA{6}(1:24); % Hd
    Pg_DA = sol_DA{13}(1:24); % Pg
    P_CHP_DA = sol_DA{8}(1:24); % PCHP


    P_DA_var = [sol_DA{3}(1:24); sol_DA{4}(1:24); Pg_DA; sol_DA{8}(1:24)];
    H_DA_var = [sol_DA{5}(1:24); sol_DA{6}(1:24); sol_DA{12}(1:24)];

    % DA performance index
    index_DA(s) = setup_cost_DA(p.v_GAS_DA,v_DA_real(times24),u_DA_real(times24),P_DA_var,H_DA_var, p.N24, p);
    index_mon_DA(s) = monetary_index_DA(p.v_GAS_DA,v_DA_real(times),u_DA_real(times),[sol_DA{13}; sol_DA{8}; sol_DA{12}], p);

    % interpolate DA solution from every hour to every 15 min
        %u_OL_rshp = [];
        Pg_DA_rshp = [];
        for jj = 1:size(u_OL,2)
          %u_OL_rshp = [u_OL_rshp repmat(u_OL(:,jj),1,Delta_DA*1/Delta_CT)];
          Pg_DA_rshp = [Pg_DA_rshp repmat(Pg_DA(jj),1,p.Delta_DA*1/p.Delta_CT)];
        end
        %u_OL_rshp = [u_OL_rshp u_OL(:,end)]; %zeros(m,Delta_DA*24/Delta_CT)];
        Pg_DA_rshp = [Pg_DA_rshp Pg_DA_rshp(end)*ones(1,p.Delta_DA/p.Delta_CT)];
        rr = rr + 96;  % trigger DA every 24h --> 24*4 (in min)

        x_OL_vec = [x_OL_vec, x_OL];
        u_OL_vec = [u_OL_vec, u_OL];
        Pg_DA_vec = [Pg_DA_vec, Pg_DA];
        P_CHP_DA_vec = [P_CHP_DA_vec, P_CHP_DA];



        % plots of the Day ahead
        set(0,'CurrentFigure',f1);

        subplot(2,4,1), cla
        hold on
        plot(1:p.N24*s+s, x_OL_vec(1,:),'b')
        grid on
        ylabel('SoC (in p.u.)')
        xlabel('Steps (in every 1h)')
        xlim([1, simEnd*p.Delta_CT])

        subplot(2,4,2), cla
        hold on
        plot(1:p.N24*s+s, x_OL_vec(2,:),'r')
        grid on
        ylabel('H (in kWh)')
        xlabel('Steps (in every 1h)')
        xlim([1, simEnd*p.Delta_CT])

        subplot(2,4,3), cla
        hold on
        plot(1:p.N24*s, u_OL_vec(1,:), 'b')
        grid on
        ylabel('Pc (in kW)')
        xlabel('Steps (in every 1h)')
        xlim([1, simEnd*p.Delta_CT])

        subplot(2,4,4), cla
        hold on
        plot(1:p.N24*s, u_OL_vec(2,:), 'b')
        grid on
        ylabel('Pd (in kW)')
        xlabel('Steps (in every 1h)')
        xlim([1, simEnd*p.Delta_CT])

        subplot(2,4,5), cla
        hold on
        plot(1:p.N24*s, u_OL_vec(3,:),'r')
        grid on
        ylabel('Hc (in kW)')
        xlabel('Steps (in every 1h)')
        xlim([1, simEnd*p.Delta_CT])

        subplot(2,4,6), cla
        hold on
        plot(1:p.N24*s, u_OL_vec(4,:),'r')
        grid on
        ylabel('Hd (in kW)')
        xlabel('Steps (in every 1h)')
        xlim([1, simEnd*p.Delta_CT])

        subplot(2,4,7), cla
        hold on
        plot(1:p.N24*s, Pg_DA_vec, 'g')
        hold on
        grid on
        ylabel('Pgrid (in kW)')
        xlabel('Steps (in every 1h)')
        xlim([1, simEnd*p.Delta_CT])

        subplot(2,4,8), cla
        hold on
        plot(1:p.N24*s, P_CHP_DA_vec, 'g')
        hold on
        grid on
        ylabel('P_{CHP} (in kW)')
        xlabel('Steps (in every 1h)')
        xlim([1, simEnd*p.Delta_CT])

        S1  = axes( 'visible', 'off', 'title', 'Day ahead simulations', 'FontSize', 16);
        drawnow

    end
     %% step 2: Continuous trading
        SoC_init_CT = x_measure(1);
        H_init_CT = x_measure(2);
        % update uncertainty profiles from the CT forecast
        times = ii:p.N_CT-1+ii;
        kk = kk+1;
        uncertainties_CT = [PL_CT(times); HL_CT(times); PV_CT(times); WT_CT(times); u_CT(times); v_CT(times)];
        % step 3: solve CT mpc problem
        mpc_CT = setup_continuous_trading(var_CT, state_CT, uncertainties_CT, Pg_DA_rshp(kk:p.N_CT-1+kk), p);
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
        [out1, out2] = dynamics(x_measure, input, p.Delta_CT); % simulate the closed-loop system
        x_measure = [out1; out2];
        x = [x, x_measure];
        u = [u, input];
        %t = [t, Delta_CT];
        Pgrid_Cls = [Pgrid_Cls, sol_CT{13}(1)];
        P_CHP_Cls = [P_CHP_Cls, sol_CT{8}(1)];
        % H_GB_Cls = [H_GB_Cls, sol_CT{12}(1)];
        Pg_DA_Cls = [Pg_DA_Cls, Pg_DA_rshp(kk)];

         % CT performance index setup_cost_CT(v_GAS_CT,v,u,P_var,H_var,Pg_DA_rshp, N_CT, Delta_CT)
        index_CT = index_CT + setup_cost_CT(p.v_GAS_CT,v_CT_real(ii),u_CT_real(ii),P_CT_var,H_CT_var, Pg_DA_rshp(kk), 1, p);
        index_mon_CT = index_mon_CT + ...
        monetary_index_CT(p.v_GAS_CT,v_CT_real(ii),u_CT_real(ii),[sol_CT{13}(1); sol_CT{8}(1); sol_CT{12}(1)], Pg_DA_rshp(kk), p);

        % plots
        set(0,'CurrentFigure',f2);

        subplot(2,4,1), cla
        hold on
        plot(x(1,:),'bo-')
        grid on
        ylabel('SoC (in p.u.)')
        xlabel('Steps (in every 15 min)')

        subplot(2,4,2), cla
        hold on
        plot(x(2,:),'ro-')
        grid on
        ylabel('H (in kWh)')
        xlabel('Steps (in every 15 min)')

        subplot(2,4,3), cla
        hold on
        plot(u(1,:),'bo-')
        grid on
        ylabel('Pc (in kW)')
        xlabel('Steps (in every 15 min)')

        subplot(2,4,4), cla
        hold on
        plot(u(2,:),'bo-')
        grid on
        ylabel('Pd (in kW)')
        xlabel('Steps (in every 15 min)')

        subplot(2,4,5), cla
        hold on
        plot(u(3,:),'ro-')
        grid on
        ylabel('Hc (in kW)')
        xlabel('Steps (in every 15 min)')

        subplot(2,4,6), cla
        hold on
        plot(u(4,:),'ro-')
        grid on
        ylabel('Hd (in kW)')
        xlabel('Steps (in every 15 min)')

        subplot(2,4,7), cla
        hold on
        plot(Pgrid_Cls,'go-')
        hold on
        plot(Pg_DA_Cls,'bo-')
        grid on
        ylabel('Pgrid_{CT}/Pgrid_{DA} (in kW)' )
        xlabel('Steps (in every 15 min)')

        subplot(2,4,8), cla
        hold on
        plot(P_CHP_Cls,'go-')
        grid on
        ylabel('P_{CHP} (in kW)')
        xlabel('Steps (in every 15 min)')

        S2  = axes( 'visible', 'off', 'title', 'Closed-loop simulations', 'FontSize', 16 );
        drawnow

 end
 toc



%% Print out the performance index
fprintf('Day ahead index is given by %3d and Continuous trading index is given by %3d\n', sum(index_DA), index_CT);
fprintf('--------------------------------------------------------------------\n');
fprintf('Day ahead monetary index is given by %3d and continuous trading monetary index is given by %3d \n', sum(index_mon_DA), index_mon_CT);

end

