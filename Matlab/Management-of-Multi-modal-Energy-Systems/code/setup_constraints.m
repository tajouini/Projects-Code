function constraints = setup_constraints(SOC_var, Heat_var, States, IC, N, p, Delta)
%% Energy Storage System (ESS)
  

    %
    % SOC_var = [I_ESS; Pc;Pd; Pload; P_WT; P_PV; P_CHP; P_EB; Pgrid];
    % Heat_var = [I_HSS; Hc; Hd; H_CHP; Hload; H_EB; H_GB];
    % State_var = [SoC; H];

    % Binary variable for charging (1) and discharging (0) (only for MIP case)
    I_ESS = SOC_var(1,:);
    % Charging power
    Pc = SOC_var(2,:);
    % Discharging power
    Pd = SOC_var(3,:);

    % Power demand
    P_load = SOC_var(4,:);
    % Wind turbine
    P_WT = SOC_var(5,:);
    % Photovoltaik
    P_PV = SOC_var(6,:);
    % Power of the combined power to heat
    P_CHP = SOC_var(7,:);
    % Heat output power of EB
	H_EB = Heat_var(6,:);
    % Power of the EB
    P_EB = SOC_var(8,:);
    % Grid power
    Pgrid = SOC_var(9,:);

    % Binary variable for HSS charging (1) and discharging (0)
    I_HSS = Heat_var(1,:);
    % Charging heat
    Hc = Heat_var(2,:);
    % Discharging heat
    Hd = Heat_var(3,:);
    % Heat output power of GB
	H_GB = Heat_var(7,:);
    % Stored heat of the combined power to heat
    H_CHP = Heat_var(4,:);
    H_load = Heat_var(5,:);

    % SOC
	SoC = States(1,:);
    % Stored hear energy
	H = States(2,:);

    %% Initial Conditions
    SoC_init = IC(1,:);
    H_init = IC(2,:);

%% Heat Storage System (HSS)
   



   %% Set up constraints
    constraints = [];
   

    constraints = [constraints, p.SocMin <= SoC(:) <= p.SocMax];
    constraints = [constraints, 0 <= Pc(:) <= p.Pmax*I_ESS(:)];
    constraints = [constraints, 0 <= Pd(:) <= p.Pmax*(1-I_ESS(:))];
    constraints = [constraints, p.Pmin_CHP <= P_CHP(:) <= p.Pmax_CHP];

    % % SOC dynamics
  
    for k=1:N
        constraints = [constraints, SoC(k+1) == p.Ad*SoC(k) + p.Bd*Delta*[Pc(k); Pd(k)]];
    end


    % HSS dynamics
   
    for k=1:N
        constraints = [constraints, H(k+1) == p.Ad_HSS*H(k) + p.Bd_HSS*Delta*[Hc(k); Hd(k)]];
    end

%% Set up HSS constraints
   

    constraints = [constraints, p.HSS_Min <= H(:) <= p.HSS_Max];
    constraints = [constraints, 0 <= Hc(:) <= p.Hmax*I_HSS(:)];
    constraints = [constraints, 0 <= Hd(:) <= p.Hmax*(1-I_HSS(:))];



   %% Set up GB constraints


    constraints = [constraints, p.H_GBMin <= H_GB(:) <= p.H_GBMax];

    %% Set up EB constraint

    constraints = [constraints, p.H_EBMin <= H_EB(:) <= p.H_EBMax];

    %% Set up network constraints

    constraints = [constraints, p.PgridMin <= Pgrid(:) <= p.PgridMax];

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
