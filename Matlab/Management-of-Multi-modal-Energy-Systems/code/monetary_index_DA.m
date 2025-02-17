function index = monetary_index_DA(v,u,var, p) % Calculate the Day Ahead Cost function

Pgrid = var(1,:);
P_CHP = var(2,:);
H_GB = var(3,:);


index = 0;

 for k=1:24

    C_e(k) = ((v(k)-u(k))/2*abs(Pgrid(k)) + (v(k)+u(k))/2*(Pgrid(k)))*p.Delta_DA;

    index = index + C_e(k);
 end

end