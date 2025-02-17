function p = load_parameters()
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    p = struct ();

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
    p.SocMin = 0.2;     % Minimal SOC % changed value
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
    p.K = 10;  % price of real time layer 
    p.alpha = 0.02;  % penalty on the input CT - real input

    p.Delta_DA = 1;
    p.Delta_CT = 1/4;



    p.Ad = 1;
    p.Bd = [p.EtaC/p.Q_ESS -1/(p.EtaD*p.Q_ESS)] ;
    % HSS dynamics   
    p.Ad_HSS = 1;
    p.Bd_HSS = [p.Eta_HC -1/p.Eta_HD];

    % initial condition 
    p.x0 = [0.3;600];









end