function [var,state] = generate_vars(N)
 % optimization variables
    SoC = sdpvar(1,N+1);
    % Binary variable for HSS charging (1) and discharging (0) 
    I_HSS = binvar(1,N);
    % Binary variable for ESS charging (1) and discharging (0) 
    I_ESS = binvar(1,N);
    % Charging heat
    Hc = sdpvar(1,N);
    % Discharging heat
    Hd = sdpvar(1,N);
    % Stored hear energy
	H = sdpvar(1,N+1);
    % Heat output power of GB
	H_GB = sdpvar(1,N);
    % Power of the combined power to heat
    P_CHP = sdpvar(1,N);
    % Electric boiler
    H_EB = sdpvar(1,N);
    % Charging heat
    Pc = sdpvar(1,N);
    % Discharging heat
    Pd = sdpvar(1,N);
    % Power network
    Pgrid = sdpvar(1,N);
    % States and variables
    var = [I_HSS; I_ESS; Hc; Hd; H_GB; P_CHP; H_EB; Pc; Pd; Pgrid];
    state = [SoC; H];
     %% initial states
   
end