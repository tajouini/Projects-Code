%%

%
% This code runs the predictive multi-energy management system.
% This implementation represents the day ahead planning finding place one
% day before the delivery and closing at 12:00- The numerical values are taken form 
% the following paper: 
%
%   "Multi-stage optimal energy amangement of multi-energy microgrid in
%    deregulated electricity markets", Wang et al., Applied
%               Energy(2022). 
%
% Copyright Institute of Automatic Control, Leibniz University Hannover.

%%
    clear; 
    clc
    close all
    yalmip('clear')
    addpath('c:\gurobi1002\win64\matlab')
%% Upper level MPC: Day ahead parameters
    Delta_t = 1;                        % sampling time in [h]
    N = 48;                             % discrete-time prediction horizon
    N24 = 24;
%% Read measurement data
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
    % Selling price (see paper)
    u_DA = 0.9*v_DA;                                         
    % purchasing price of natural gas (see paper)
    v_GAS_DA = 0.055;   
    % % Real purchasing price 
    load('C:\Users\Jouini\Documents\Research\IFES-Multi-modale-Energy-systeme\multi-energy-mangement-mpc\data\8760h-01p\v_DA_real.mat');
    u_DA_real = 0.91*v_DA_real;

    
    % Power measurements of WT generation 
    P_WT = WT_DA(1:2*N24);     % in kW
    % Power measurements of PV generation 
    P_PV = PV_DA(1:2*N24);
    % Power demand 
    P_load = PL_DA(1:2*N24);    % in kW
    % Heat demand 
    H_load = HL_DA(1:2*N24);
    % Price measurements 
    v_data = v_DA(1:2*N24);                                  % purchasing price of electricity
    u_data = 0.9*v_data;                                         % Selling price (see paper)
    v_GAS = 0.055;                                               % purchasing price of natural gas (see paper)
        
%% Energy Storage System (ESS)
% Set up ESS parameters 
    Pmax = 500;                        % Maximum charging/discharging power
    % Capacity of ESS
    Q_ESS = 2000;                      % in [kWh]
    % Charging efficiency
    EtaC = 0.95;                       % in p.u.
    % Discharging efficiency
    EtaD = 0.95;                       % in p.u.
    % Depreciation coefficient         % in $/kWh
    Rho_ESS = 0.01;
    % Characteristic of the heat-to-electric ratio
    b = 1.2;

    % Binary variable for charging (1) and discharging (0) (only for MIP case)
    I_ESS = binvar(1,N);
    % Charging power
    Pc = sdpvar(1,N);
    % Discharging power
    Pd = sdpvar(1,N);
    % SOC
	SoC = sdpvar(1,N+1);
    % initial SOC
    SoC_init = sdpvar(1,1);
    

%% Heat Storage System (HSS)
    % Set up ESS parameters 
    % Charging efficiency
    Eta_HC = 0.9;                       % in p.u.
    % Discharging efficiency
    Eta_HD = 0.9;                       % in p.u.
    % Depreciation coefficient          % in $/kWh
    Rho_HSS = 0.01;
    % Gas to electricity efficiency of CHP
    Eta_CHP = 0.95;                     % in p.u.
    % Binary variable for HSS charging (1) and discharging (0) 
    I_HSS = binvar(1,N);
    % Charging heat
    Hc = sdpvar(1,N);
    % Discharging heat
    Hd = sdpvar(1,N);
    % Stored hear energy
	H = sdpvar(1,N+1);
    % Initial stored heat energy
    H_init = sdpvar(1,1);
    % Power of the combined power to heat
    P_CHP = sdpvar(1,N);
    % Stored heat of the combined power to heat
    H_CHP = b*P_CHP; 
    
%% Electric Boiler System (EB) - Electricity to heat
    % Heat output power of EB
	H_EB = sdpvar(1,N);
    % Electricity to heat efficiency of the EB
    Eta_EB = 0.95;   
    % Power of the EB 
    P_EB = H_EB/Eta_EB;

    
%% Gas Boiler System (GB) - Gas to heat
    % Electricity to heat efficiency of the GB
    Eta_GB = 0.8;
    % Heat output power of GB
	H_GB = sdpvar(1,N);
    
%% Network
    % Grid power
    Pgrid   = sdpvar(1,N);

%% Set up ESS constraints
    constraints = [];
    % Maximal SOC 
    SocMax = 0.8;
    % Minimal SOC 
    SocMin = 0.2;
    % Maximal combined power to heat (in kW)
    Pmax_CHP = 1200;
    % minimal combined power to heat (in kW)
    Pmin_CHP = 0;

    constraints = [constraints, SocMin <= SoC(:) <= SocMax];
    constraints = [constraints, 0 <= Pc(:) <= Pmax*I_ESS(:)];
    constraints = [constraints, 0 <= Pd(:) <= Pmax*(1-I_ESS(:))];
    constraints = [constraints, Pmin_CHP <= P_CHP(:) <= Pmax_CHP];

    % SOC dynamics   
    Ad = 1;
    Bd = [EtaC/Q_ESS -1/(EtaD*Q_ESS)]*Delta_t;
    
    for k=1:N
        constraints = [constraints, SoC(k+1) == Ad*SoC(k) + Bd*[Pc(k); Pd(k)]];
    end

%% Set up HSS constraints
    % Maximal stored heat (in kWh)
    HSS_Max = 3000;
    % Minimal stored heat (in kWh)
    HSS_Min = 100;
    % Maximum chargin/discharding power
    Hmax = 750;

    constraints = [constraints, HSS_Min <= H(:) <= HSS_Max];
    constraints = [constraints, 0 <= Hc(:) <= Hmax*I_HSS(:)];
    constraints = [constraints, 0 <= Hd(:) <= Hmax*(1-I_HSS(:))];
   
    % HSS dynamics   
    Ad_HSS = 1;
    Bd_HSS = [Eta_HC -1/Eta_HD]*Delta_t;
    
    for k=1:N
        constraints = [constraints, H(k+1) == Ad_HSS*H(k) + Bd_HSS*[Hc(k); Hd(k)]];
    end
    
   %% Set up GB constraints
    % Maximal heat power limit of the GB
    H_GBMax = 500;
    % Minimal heat power limit of the GB
    H_GBMin = 0; % I changed this lower bound
    

    constraints = [constraints, H_GBMin <= H_GB(:) <= H_GBMax];
   
    %% Set up EB constraints
    % Maximal heat power limit of the EB
    H_EBMax = 500;
    % Minimal heat power limit of the EB
    H_EBMin = 0;
    
    constraints = [constraints, H_EBMin <= H_EB(:) <= H_EBMax];

    %% Set up network constraints
    % Maximal grid power
    PgridMax = 2500;
    % Minimal grid power
    PgridMin = -2500;
    
    constraints = [constraints, PgridMin <= Pgrid(:) <= PgridMax];
    
 
    %% Heat and Power balance 
    constraints = [constraints, Pgrid(:)+P_PV(:)+P_WT(:)+Pd(:)+P_CHP(:) == Pc(:)+P_EB(:)+P_load(:)];
    % Heat balance 
    
    constraints = [constraints, H_CHP(:) + H_GB(:) + H_EB(:) + Hd(:) == Hc(:)+ H_load(:)];
    
    % Purchasing price
    v = v_data(1:N); 
    % Selling price
    u = u_data(1:N); 
    
    C_e = sdpvar(1,N);
    C_ESS = sdpvar(1,N);
    C_HSS = sdpvar(1,N);
    C_GAS = sdpvar(1,N);


%% Day Ahead Cost function 
    cost = 0;
    for k=1:N
    C_e(k) = ((v(k)-u(k))/2*abs(Pgrid(k)) + (v(k)+u(k))/2*(Pgrid(k)))*Delta_t;
    C_ESS(k) =  Rho_ESS*(Pd(k)+Pc(k))*Delta_t;
    C_HSS(k) =  Rho_HSS*(Hd(k)+Hc(k))*Delta_t;
    C_GAS(k) = v_GAS * (P_CHP(k)/Eta_CHP+H_GB(k)/Eta_GB);

    cost = cost + C_e(k) + C_ESS(k) + C_HSS(k) + C_GAS(k);
    end
    
    constraints = [constraints, SoC(1) == SoC_init];
    constraints = [constraints, H(1) == H_init];
   
    % solver settings
    ops = sdpsettings('solver','gurobi','verbose',5);
    % Initial conditions/predictions
    init = {SoC_init, H_init};
    % Output of 'mpc' after solving the optimization problem
    out = {SoC,Pc,Pd,P_CHP,P_EB,I_ESS, H, Hc, Hd, I_HSS, H_EB, H_GB, Pgrid, cost};
    % Generate yalmip object for optimization in closed-loop
    mpc = optimizer(constraints,cost,ops,init,out);    

    % initial soc is 0.3 and initial storage is 600 (see paper)
    [sol,flag] = mpc({0.3,600});

    SoC_Sol = sol{1};
    Pc_Sol = sol{2};
    Pd_Sol = sol{3};
    P_CHP_Sol = sol{4};
    P_EB_Sol = sol{5};
    I_ESS_Sol = sol{6};
    H_Sol = sol{7};
    Hc_Sol = sol{8};
    Hd_Sol = sol{9};
    I_HSS_Sol = sol{10};
    H_EB_Sol = sol{11};
    H_GB_Sol = sol{12};
    Pgrid_Sol = sol{13};
    cost_opt = sol{14};
    
    cost24 = 0;
    for k=1:24
    C_e(k) = ((v(k)-u(k))/2*abs(Pgrid_Sol(k)) + (v(k)+u(k))/2*(Pgrid_Sol(k)))*Delta_t;
    C_ESS(k) =  Rho_ESS*(Pd_Sol(k)+Pc_Sol(k))*Delta_t;
    C_HSS(k) =  Rho_HSS*(Hd_Sol(k)+Hc_Sol(k))*Delta_t;
    C_GAS(k) = v_GAS * (P_CHP_Sol(k)/Eta_CHP+H_GB_Sol(k)/Eta_GB);

    cost24 = cost24 + C_e(k) + C_ESS(k) + C_HSS(k) + C_GAS(k);
    end

    P_DA_var = [Pc_Sol(1:24); Pd_Sol(1:24); Pgrid_Sol(1:24); P_CHP_Sol(1:24)];
    H_DA_var = [Hc_Sol(1:24); Hd_Sol(1:24); H_GB_Sol(1:24)];
    index_DA = setup_cost_DA(v_GAS_DA,v_DA_real(1:N24),u_DA_real(1:N24),P_DA_var,H_DA_var, N24, Delta_t);

%% Plots
% Energy system
figure(1)

ax1 = subplot(3,2,1);
plot(Pc_Sol,'k--')
hold on
plot(Pd_Sol,'g--')
grid on
legend('Pc (in kW)','Pd (in kW)')
xlim([1,24])

ax2 = subplot(3,2,2);
plot(P_CHP_Sol, 'r')
hold on
grid on
legend('P_{CHP} (in kW)')
xlim([1,24])

ax3 = subplot(3,2,3);
plot(Pgrid_Sol,'b')
grid on
legend('Pgrid (in kW)')
xlim([1,24])

ax4 = subplot(3,2,4);
plot(SoC_Sol, 'b')
grid on
legend('Soc (in p.u.)');
xlim([1,24])

ax5 = subplot(3,2,5);
plot(I_ESS_Sol)
grid on
legend('I_{ESS}');
xlim([1,24])

ax5 = subplot(3,2,6);
plot(cost24*ones(1,24))
grid on
legend('Cost of 24h');
xlim([1,24])

sgtitle('48h simulation of ESS')

% Heat system
figure(2)

ax1 = subplot(3,2,1);
plot(Hc_Sol,'k--')
grid on
hold on 
plot(Hd_Sol, 'g')
legend('H_c (in kW)')
xlim([1,24])

% ax2 = subplot(3,2,2);
% plot(Hd_Sol, 'g')
% grid on
% legend('Hd in kWh');
% xlim([1,24])

ax3 = subplot(3,2,2);
plot(H_EB_Sol, 'r')
grid on
legend('H_{EB} (in kW)');
xlim([1,24])

ax4 = subplot(3,2,3);
plot(H_GB_Sol, 'b')
grid on
legend('H_{GB} (in kW)');
xlim([1,24])


ax5 = subplot(3,2,4);
plot(I_HSS_Sol)
grid on
legend('I_{HSS}');
xlim([1,24])

ax6 = subplot(3,2,5);
plot(b*P_CHP_Sol)
grid on
legend('H_{CHP}');
xlim([1,24])


ax7 = subplot(3,2,6);
plot(H_Sol, 'g')
grid on
legend('H in kWh');
xlim([1,24])


sgtitle('48h simulation of gas and HSS')

% load and generation profiles
% load and generation profiles
figure(3)

ax1 = subplot(4,1,1);
plot(v)
hold on 
plot(u)
grid on
title('Prices')
xlabel('Time in hours')
ylabel('Price in dollars')
xlim([1,48])

ax2 = subplot(4,1,2);
plot(P_load);
grid on
hold on
plot(H_load);
title('Load profile')
xlabel('Time in hours')
ylabel('in kW')
xlim([1,48])

ax3 = subplot(4,1,3);
plot(P_PV);
grid on
title('PV generation')
xlabel('Time in hours')
ylabel('Power in kW')
xlim([1,48])


ax4 = subplot(4,1,4);
plot(P_WT);
grid on
title('Wind Turbine generation')
xlabel('Time in hours')
ylabel('Power in kW')
xlim([1,48])

sgtitle('48h Uncertainty profiles')

% power and heat balance
figure(4)

% check power balance Pgrid+P_PV+P_WT+Pd+P_CHP = Pc+P_EB+P_load
ax41 = subplot(2,1,1);
plot(Pgrid_Sol+P_PV+P_WT+Pd_Sol+P_CHP_Sol)
hold on 
plot(Pc_Sol+P_EB_Sol+P_load)
grid on
xlabel('Time in hours')
ylabel('Power in kW')
legend('Power Generation','Power Demand')

% check heat balance H_CHP + H_GB + H_EB + Hd = Hc+ H_load
x42 = subplot(2,1,2);
plot(b*P_CHP_Sol+ H_GB_Sol + H_EB_Sol + Hd_Sol)
hold on 
plot(Hc_Sol+ H_load)
grid on
xlabel('Time in hours')
ylabel('Power in kW')
legend('Power heat Generation','Power heat Demand')


sgtitle('Heat and power balance')

