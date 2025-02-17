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
 
function mpc_CT = setup_continuous_trading_v2(var,state, uncertainties_CT, Pg_DA_rshp,v_GAS_CT, N_CT, Delta_CT)
 
   
    
    %% data measurements 
    Pload = uncertainties_CT(1,:);
    Hload = uncertainties_CT(2,:);  
    P_PV = uncertainties_CT(3,:);
    P_WT = uncertainties_CT(4,:);
    u = uncertainties_CT(5,:); 
    v = uncertainties_CT(6,:); 
                                          
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
    
    Pgrid = P_v + P_u + Pg_DA_rshp;
   
    % setup the constraints
    SOC_var = [Pc;Pd; Pload; P_WT; P_PV; P_CHP; P_EB; P_v; P_u];
    Heat_var = [Hc; Hd; H_CHP; Hload; H_EB; H_GB];
    State_var = [SoC; H];
    IC = [SoC_init; H_init];
    constraints_CT = setup_constraints_v2(SOC_var, Heat_var, State_var,IC, N_CT, Delta_CT,"CT", Pg_DA_rshp);
    

    % setup the cost of DA
    bin_var = [beta_ESS; gamma_HSS];
    P_var = [Pc;Pd;P_v;P_u;P_CHP];
    H_var = [Hc; Hd; H_GB];
    cost_CT = setup_cost_CT_v2(v_GAS_CT,v,u,bin_var,P_var,H_var, N_CT, Delta_CT); 

    % solver settings
    ops = sdpsettings('solver','gurobi','verbose',5);
    % Initial conditions/predictions
    init = {SoC_init, H_init};
    % Output of 'mpc' after solving the optimization problem
    out = {SoC,H,Pc,Pd,Hc,Hd,P_EB,P_CHP,beta_ESS,gamma_HSS,H_EB,H_GB,Pgrid,cost_CT};
    % Generate yalmip object for optimization in closed-loop
    mpc_CT = optimizer(constraints_CT,cost_CT,ops,init,out);       

end    
