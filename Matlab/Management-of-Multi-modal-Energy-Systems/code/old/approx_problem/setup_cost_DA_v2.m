function cost_DA = setup_cost_DA_v2(v_GAS_DA,v,u,bin_var, P_var,H_var, N_DA, Delta_DA) % Calculate the Day Ahead Cost function

% P_var = [Pc;Pd;Pgrid;P_CHP];
% H_var = [Hc; Hd; H_GB];
beta_ESS = bin_var(1,:);
gamma_HSS = bin_var(2,:);
Pc= P_var(1,:);
Pd = P_var(2,:);
P_v = P_var(3,:);
P_u = P_var(4,:);
P_CHP = P_var(5,:);
Hc = H_var(1,:);
Hd = H_var(2,:);
H_GB = H_var(3,:);
% Depreciation coefficient         % in $/kWh
Rho_ESS = 0.01;
% Depreciation coefficient          % in $/kWh
Rho_HSS = 0.01;
% purchasing price of gas
v_GAS = v_GAS_DA;   
% Gas to electricity efficiency of CHP
Eta_CHP = 0.95; 
% Electricity to heat efficiency of the GB
Eta_GB = 0.8;
% epsilon 
eps = 1e-4;

cost_DA = 0;

 for k=1:N_DA
    C_e_new(k) = (v(k)*P_v(k)+u(k)*P_u(k))*Delta_DA;
    C_ESS(k) =  Rho_ESS*(Pd(k)+Pc(k))*Delta_DA;
    C_HSS(k) =  Rho_HSS*(Hd(k)+Hc(k))*Delta_DA;
    C_GAS(k) = v_GAS*(P_CHP(k)/Eta_CHP + H_GB(k)/Eta_GB);
    C_bin_ESS(k) = (1/eps*beta_ESS(k)+eps*(Pc(k)-Pd(k)))^2*Delta_DA;
    C_bin_HSS(k) = (1/eps*gamma_HSS(k)+eps*(Hc(k)-Hd(k)))^2*Delta_DA;
    cost_DA = cost_DA + C_e_new(k) + C_ESS(k) + C_HSS(k) + C_GAS(k) + C_bin_ESS(k) + C_bin_HSS(k);
 end


end