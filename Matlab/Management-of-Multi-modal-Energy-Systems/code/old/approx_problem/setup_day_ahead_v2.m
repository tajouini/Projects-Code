%%

%
% This code sets up the day ahead planning finding place one
% day before the delivery and closing at 12:00- The numerical values are taken form 
% the following paper: 
%
%   "Multi-stage optimal energy amangement of multi-energy microgrid in
%    deregulated electricity markets", Wang et al., Applied
%               Energy(2022). 
%
% Copyright Institute of Automatic Control, Leibniz University Hannover.



function mpc_day_ahead = setup_day_ahead_v2(var,state, uncertainties_DA, v_GAS_DA,N_DA, Delta_DA)

   
    %% parameters
    % Day ahead parameters
    N = N_DA;                             % discrete-time prediction horizon
   
    % data measurements 
    Pload = uncertainties_DA(1,:);
    Hload = uncertainties_DA(2,:);  
    P_PV = uncertainties_DA(3,:);
    P_WT = uncertainties_DA(4,:);
    u = uncertainties_DA(5,:); 
    v = uncertainties_DA(6,:); 
                                          
   
    %% initial states
    SoC_init = sdpvar(1,1);
    H_init = sdpvar(1,1);

    %% optimization variables
    %  % States and variables
    % var = [gamma_HSS; beta_ESS; Hc; Hd; H_GB; P_CHP; H_EB; Pc; Pd; Pgrid];
    % state = [SoC; H];
    SoC = state(1,:);
    % Stored hear energy
	H = state(2,:);
    % Binary variable for HSS charging (1) and discharging (0) 
    gamma_HSS = var(1,:);
    % Binary variable for ESS charging (1) and discharging (0) 
    beta_ESS = var(2,:);
    % Charging heat
    Hc = var(3,:);
    % Discharging heat
    Hd = var(4,:);
   
    % Heat output power of GB
	H_GB = var(5,:);
    % Power of the combined power to heat
    P_CHP = var(6,:);
    % Characteristic of the heat-to-electric ratio
    b = 1.2;
    H_CHP = b*P_CHP; 
    % Electric boiler
    H_EB = var(7,:);
    % Electricity to heat efficiency of the EB
    Eta_EB = 0.95;   
    % Power of the EB 
    P_EB = H_EB/Eta_EB;
    % Charging heat
    Pc = var(8,:);
    % Discharging heat
    Pd = var(9,:);
    % Power network
    P_v = var(10,:);
    P_u = var(11,:);

    Pgrid = P_v + P_u;
    
   
    % setup the constraints
    SOC_var = [Pc;Pd; Pload; P_WT; P_PV; P_CHP; P_EB; P_v; P_u];
    Heat_var = [Hc; Hd; H_CHP; Hload; H_EB; H_GB];
    State_var = [SoC; H];
    IC = [SoC_init; H_init];
    constraints_DA = setup_constraints_v2(SOC_var, Heat_var, State_var,IC, N_DA, Delta_DA,"DA",[]);
   

    % setup the cost of DA
    bin_var = [beta_ESS; gamma_HSS];
    P_var = [Pc;Pd;P_v;P_u;P_CHP];
    H_var = [Hc; Hd; H_GB];
    cost_DA = setup_cost_DA_v2(v_GAS_DA,v,u,bin_var,P_var,H_var, N_DA, Delta_DA); 

    % solver settings
    ops = sdpsettings('solver','gurobi','verbose',5);
    % Initial conditions/predictions
    init = {SoC_init, H_init};
    % Output of 'mpc' after solving the optimization problem
    out = {SoC,H,Pc,Pd,Hc,Hd,P_EB,P_CHP,beta_ESS,gamma_HSS,H_EB,H_GB,Pgrid,cost_DA};
    % Generate yalmip object for optimization in closed-loop
    mpc_day_ahead = optimizer(constraints_DA,cost_DA,ops,init,out);       
end    
