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



function mpc_day_ahead = setup_day_ahead(var,state, uncertainties_DA, N_DA, p)

   
    %% parameters
    % Day ahead parameters
  
   
    % data measurements 
    Pload = uncertainties_DA(1,:);
    Hload = uncertainties_DA(2,:);  
    P_PV = uncertainties_DA(3,:);
    P_WT = uncertainties_DA(4,:);
    u = uncertainties_DA(5,:); 
    v = uncertainties_DA(6,:); 
                                          
   %% optimization variables
    %  % States and variables
    % var = [I_HSS; I_ESS; Hc; Hd; H_GB; P_CHP; H_EB; Pc; Pd; Pgrid];
    % state = [SoC; H];
    SoC = state(1,:);
    % Stored hear energy
	H = state(2,:);
    % Binary variable for HSS charging (1) and discharging (0) 
    I_HSS = var(1,:);
    % Binary variable for ESS charging (1) and discharging (0) 
    I_ESS = var(2,:);
    % Charging heat
    Hc = var(3,:);
    % Discharging heat
    Hd = var(4,:);
   
    % Heat output power of GB
	H_GB = var(5,:);
    % Power of the combined power to heat
    P_CHP = var(6,:);
    H_CHP = p.b*P_CHP; 
    % Electric boiler
    H_EB = var(7,:);
    % Power of the EB 
    P_EB = H_EB/p.Eta_EB;
    % Charging heat
    Pc = var(8,:);
    % Discharging heat
    Pd = var(9,:);
    % Power network
    Pgrid = var(10,:);
    
    SoC_init_DA = sdpvar(1,1);
    H_init_DA = sdpvar(1,1);
    IC = [SoC_init_DA ; H_init_DA];
    % setup the constraints
    SOC_var = [I_ESS; Pc;Pd; Pload; P_WT; P_PV; P_CHP; P_EB; Pgrid];
    Heat_var = [I_HSS; Hc; Hd; H_CHP; Hload; H_EB; H_GB];
    State_var = [SoC; H];
    constraints_DA = setup_constraints(SOC_var, Heat_var, State_var,IC, N_DA, p, p.Delta_DA);
   

    % setup the cost of DA
    P_var = [Pc;Pd;Pgrid;P_CHP];
    H_var = [Hc; Hd; H_GB];
    % setup_cost_DA(v_GAS_DA,v,u,P_var,H_var, N_DA, Delta_DA)
    cost_DA = setup_cost_DA(v,u,P_var,H_var, N_DA, p); 

    % solver settings
    if exist ('OCTAVE_HOME','builtin') == 5
      ops = sdpsettings('solver','glpk');
    else
      ops = sdpsettings('solver','gurobi','verbose',5);
    end    % Initial conditions/predictions
    % Initial conditions/predictions
    init = {SoC_init_DA, H_init_DA};
    % Output of 'mpc' after solving the optimization problem
    out = {SoC,H,Pc,Pd,Hc,Hd,P_EB,P_CHP,I_ESS,I_HSS,H_EB,H_GB,Pgrid,cost_DA};
    % Generate yalmip object for optimization in closed-loop
    mpc_day_ahead = optimizer(constraints_DA,cost_DA,ops,init,out);       
end    
