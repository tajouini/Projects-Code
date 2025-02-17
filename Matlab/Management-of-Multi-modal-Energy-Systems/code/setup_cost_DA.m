function cost_DA = setup_cost_DA(v,u,P_var,H_var, N_DA, p) % Calculate the Day Ahead Cost function

% P_var = [Pc;Pd;Pgrid;P_CHP];
% H_var = [Hc; Hd; H_GB];

Pc= P_var(1,:);
Pd = P_var(2,:);
Pgrid = P_var(3,:);
P_CHP = P_var(4,:);
Hc = H_var(1,:);
Hd = H_var(2,:);
H_GB = H_var(3,:);


cost_DA = 0;

 for k=1:N_DA

    C_e(k) = ((v(k)-u(k))/2*abs(Pgrid(k)) + (v(k)+u(k))/2*(Pgrid(k)))*p.Delta_DA;
    C_ESS(k) =  p.Rho_ESS*(Pd(k)+Pc(k))*p.Delta_DA;
    C_HSS(k) =  p.Rho_HSS*(Hd(k)+Hc(k))*p.Delta_DA;
    C_GAS(k) = p.v_GAS_DA*(P_CHP(k)/p.Eta_CHP + H_GB(k)/p.Eta_GB)*p.Delta_DA;

    cost_DA = cost_DA + C_e(k) + C_ESS(k) + C_HSS(k) + C_GAS(k);
 end

end