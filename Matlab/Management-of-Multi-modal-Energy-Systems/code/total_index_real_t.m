function [C_tot, C_input,  C_e] = total_index_real_t(u_real, input_CT, P_CHP, H_GB, Pgrid_real,Pg_CT, p) 
% Optimization variables
Pc = u_real(1);
Pd = u_real(2);
Hc = u_real(3);
Hd = u_real(4);

C_e = abs(Pgrid_real - Pg_CT)*p.Delta_CT*p.K;
C_ESS =  p.Rho_ESS*(Pd + Pc)*p.Delta_CT;
C_HSS =  p.Rho_HSS*(Hd + Hc)*p.Delta_CT;
C_GAS = p.v_GAS_CT*(P_CHP/p.Eta_CHP + H_GB/p.Eta_GB)*p.Delta_CT;
C_input = p.alpha* norm(input_CT-u_real,2);

C_tot = C_e + C_ESS + C_HSS + C_GAS; 
%C_index2 = ;

end