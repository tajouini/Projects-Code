function cost_real_t = setup_cost_real_t(v,u,P_var,H_var,Pg_CT, input_CT,N_real_t,p) % Calculate the real time Cost function

% Optimization variables
Pc= P_var(1,:);
Pd = P_var(2,:);
Pgrid_real = P_var(3,:);
P_CHP = P_var(4,:);
Hc = H_var(1,:);
Hd = H_var(2,:);
H_GB = H_var(3,:);
input = [Pc,Pd,Hc,Hd];

cost_real_t = 0;

% cost 1

    C_e = abs(Pgrid_real - Pg_CT)*p.Delta_CT*p.K;
    C_ESS =  p.Rho_ESS*(Pd + Pc)*p.Delta_CT;
    C_HSS =  p.Rho_HSS*(Hd + Hc)*p.Delta_CT;
    C_GAS = p.v_GAS_CT*(P_CHP/p.Eta_CHP + H_GB/p.Eta_GB)*p.Delta_CT;


    cost_real_t = cost_real_t + C_e + C_ESS + C_HSS + C_GAS;

    C_input = p.alpha* norm(input_CT-input,2);

    cost_real_t = cost_real_t+ C_input;

 % cost 2 cost_real = cost_CT


end