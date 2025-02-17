function index = monetary_index_CT(v,u,var,Pg_DA_rshp, N_CT,p) % Calculate the Day Ahead Cost function

% Optimization variables

Pgrid = var(1,:);
%P_CHP = var(2,1);
%H_GB = var(3,1);

index = 0;

for k=1:N_CT
C_e = ((v(k)-u(k))/2*abs(Pgrid(k)-Pg_DA_rshp(k)) + (v(k)+u(k))/2*(Pgrid(k)-Pg_DA_rshp(k)))*p.Delta_CT;
%C_GAS = p.v_GAS_CT*(P_CHP/p.Eta_CHP + H_GB/p.Eta_GB)*p.Delta_CT;
index =  index + C_e; %+ C_GAS;
end

end