function constraints = setup_constraints_v2(SOC_var, Heat_var, States, IC, N, Delta_t, letter, Pg_DA_rshp)
%% Energy Storage System (ESS)
    % Set up ESS parameters 
    Pmax = 500;                        % Maximum charging/discharging power
    % Capacity of ESS
    Q_ESS = 2000;                      % in [kWh]
    % Charging efficiency
    EtaC = 0.95;                       % in p.u.
    % Discharging efficiency
    EtaD = 0.95;                       % in p.u.

    % 
    % SOC_var = [Pc;Pd; Pload; P_WT; P_PV; P_CHP; P_EB; P_v; P_u];
    % Heat_var = [Hc; Hd; H_CHP; Hload; H_EB; H_GB];
    % State_var = [SoC; H];

    % Binary variable for charging (1) and discharging (0) (only for MIP case)
    %beta_ESS = SOC_var(1,:);
    % Charging power
    Pc = SOC_var(1,:);
    % Discharging power
    Pd = SOC_var(2,:);
    
    % Power demand
    P_load = SOC_var(3,:);
    % Wind turbine
    P_WT = SOC_var(4,:);
    % Photovoltaik
    P_PV = SOC_var(5,:);
    % Power of the combined power to heat
    P_CHP = SOC_var(6,:);
    % Power of the EB 
    P_EB = SOC_var(7,:);
    % Grid power
    P_v = SOC_var(8,:);
    P_u = SOC_var(9,:);
    if letter=="DA"
        Pgrid = P_v + P_u;
    else 
        if letter=="CT"
        Pgrid = P_v + P_u + Pg_DA_rshp;
        end
    end
    % Binary variable for HSS charging (1) and discharging (0) 
    %gamma_HSS = Heat_var(1,:);
    % Charging heat
    Hc = Heat_var(1,:);
    % Discharging heat
    Hd = Heat_var(2,:);
    % Stored heat of the combined power to heat
    H_CHP = Heat_var(3,:);
    H_load = Heat_var(4,:);
    % Heat output power of EB
	H_EB = Heat_var(5,:);
    % Heat output power of GB
	H_GB = Heat_var(6,:);
   
    % SOC
	SoC = States(1,:);
    % Stored hear energy
	H = States(2,:);

    %% Initial Conditions
    SoC_init = IC(1,:);
    H_init = IC(2,:);

%% Heat Storage System (HSS)
    % Charging efficiency
    Eta_HC = 0.9;                       % in p.u.
    % Discharging efficiency
    Eta_HD = 0.9;                       % in p.u.
   
    
    
   
   %% Set up constraints
    constraints = [];
  % for reformulated cost
    constraints = [constraints, 0 <= P_v(:)];
    constraints = [constraints, P_u(:) <= 0];

    % Maximal SOC 
    SocMax = 0.8;
    % Minimal SOC 
    SocMin = 0.2;
    % Maximal combined power to heat (in kW)
    Pmax_CHP = 1200;
    % minimal combined power to heat (in kW)
    Pmin_CHP = 0;

    constraints = [constraints, SocMin <= SoC(:) <= SocMax];
    constraints = [constraints, 0 <= Pc(:) <= Pmax];
    constraints = [constraints, 0 <= Pd(:) <= Pmax];
    constraints = [constraints, Pmin_CHP <= P_CHP(:) <= Pmax_CHP];

    % % SOC dynamics   
    Ad = 1;
    Bd = [EtaC/Q_ESS -1/(EtaD*Q_ESS)]*Delta_t;

    for k=1:N
        constraints = [constraints, SoC(k+1) == Ad*SoC(k) + Bd*[Pc(k); Pd(k)]];
    end


    % HSS dynamics   
    Ad_HSS = 1;
    Bd_HSS = [Eta_HC -1/Eta_HD]*Delta_t;
    for k=1:N
        constraints = [constraints, H(k+1) == Ad_HSS*H(k) + Bd_HSS*[Hc(k); Hd(k)]];
    end
    
%% Set up HSS constraints
    % Maximal stored heat (in kWh)
    HSS_Max = 3000;
    % Minimal stored heat (in kWh)
    HSS_Min = 100;
    % Maximum chargin/discharding power
    Hmax = 750;

    constraints = [constraints, HSS_Min <= H(:) <= HSS_Max];
    constraints = [constraints, 0 <= Hc(:) <= Hmax];
    constraints = [constraints, 0 <= Hd(:) <= Hmax];
   
   
    
   %% Set up GB constraints
    % Maximal heat power limit of the GB
    H_GBMax = 500;
    % Minimal heat power limit of the GB
    H_GBMin = 0;
    

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
    
    % Power balance 
    constraints = [constraints, Pgrid(:) + P_PV(:) + P_WT(:) + Pd(:) + P_CHP(:) == ...
        Pc(:) + P_EB(:) + P_load(:)];
    % Heat balance 
    constraints = [constraints, H_CHP(:) + H_GB(:) + H_EB(:) + Hd(:) ==...
        Hc(:) + H_load(:)];


    % initialization
    constraints = [constraints, SoC(1) == SoC_init];
    constraints = [constraints, H(1) == H_init];

 end