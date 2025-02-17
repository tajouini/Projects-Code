function cost_DA = monetary_index_DA_v2(v_GAS_DA,v,u,P_var,H_var, Delta_DA) % Calculate the Day Ahead Cost function

P_v = P_var(3,:);
P_u = P_var(4,:);
P_CHP = P_var(5,:);
H_GB = H_var(3,:);

% purchasing price of gas
v_GAS = v_GAS_DA;   
% Gas to electricity efficiency of CHP
Eta_CHP = 0.95; 
% Electricity to heat efficiency of the GB
Eta_GB = 0.8;

cost_DA = 0;

 for k=1:24
    C_e_new(k) = (v(k)*P_v(k)+u(k)*P_u(k))*Delta_DA;
    C_GAS(k) = v_GAS*(P_CHP(k)/Eta_CHP + H_GB(k)/Eta_GB);
    cost_DA = cost_DA + C_e_new(k) + C_GAS(k);
 end


end