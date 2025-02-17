function cost_CT = setup_cost_CT(v,u,P_var,H_var,Pg_DA_rshp, N_CT,p) % Calculate the Day Ahead Cost function

% Optimization variables
Pc= P_var(1,:);
Pd = P_var(2,:);
Pgrid = P_var(3,:);
P_CHP = P_var(4,:);
Hc = H_var(1,:);
Hd = H_var(2,:);
H_GB = H_var(3,:);


cost_CT = 0;

 for k=1:N_CT
   C_e = ((v(k)-u(k))/2*abs(Pgrid(k)-Pg_DA_rshp(k)) + (v(k)+u(k))/2*(Pgrid(k)-Pg_DA_rshp(k)))*p.Delta_CT;
    C_ESS =  p.Rho_ESS*(Pd(k)+Pc(k))*p.Delta_CT;
    C_HSS =  p.Rho_HSS*(Hd(k)+Hc(k))*p.Delta_CT;
    C_GAS = p.v_GAS_CT*(P_CHP(k)/p.Eta_CHP + H_GB(k)/p.Eta_GB)*p.Delta_CT;

    cost_CT = cost_CT + C_e + C_ESS + C_HSS + C_GAS;
 end

end