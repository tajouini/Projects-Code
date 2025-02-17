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
 
function mpc_CT = setup_real_t(var,state, real_time_profiles, Pg_CT,input_CT, N_real_t,p)
 
   
    
    %% data measurements 
    Pload = real_time_profiles(1,:);
    Hload = real_time_profiles(2,:);  
    P_PV = real_time_profiles(3,:);
    P_WT = real_time_profiles(4,:);
    u = real_time_profiles(5,:); 
    v = real_time_profiles(6,:); 
                                          
    %% initial states
    SoC_init = sdpvar(1,1);
    H_init = sdpvar(1,1);

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
    
   
    % setup the constraints
    SOC_var = [I_ESS; Pc;Pd; Pload; P_WT; P_PV; P_CHP; P_EB; Pgrid];
    Heat_var = [I_HSS; Hc; Hd; H_CHP; Hload; H_EB; H_GB];
    State_var = [SoC; H];
    IC = [SoC_init; H_init];
    constraints_real_t = setup_constraints(SOC_var, Heat_var, State_var,IC, ...
        N_real_t, p, p.Delta_CT);
   

    % setup the cost of DA
    P_var = [Pc;Pd;Pgrid;P_CHP];
    H_var = [Hc; Hd; H_GB];
    % (K,P_var,H_var,Pg_CT,p)
    cost_CT = setup_cost_real_t(v,u,P_var,H_var,Pg_CT,input_CT,N_real_t, p); 

    % solver settings
    if exist ('OCTAVE_HOME','builtin') == 5
      ops = sdpsettings('solver','glpk');
    else
      ops = sdpsettings('solver','gurobi','verbose',5);
    end    % Initial conditions/predictions
    init = {SoC_init, H_init};
    % Output of 'mpc' after solving the optimization problem
    out = {SoC,H,Pc,Pd,Hc,Hd,P_EB,P_CHP,I_ESS,I_HSS,H_EB,H_GB,Pgrid,H_CHP};
    % Generate yalmip object for optimization in closed-loop
    mpc_CT = optimizer(constraints_real_t,cost_CT,ops,init,out);       

end    
