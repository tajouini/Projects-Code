function [SoC_next, H_next] = dynamics(IC, input, p)
SoC_0 = IC(1);
H_SS_0 = IC(2);
Pc = input(1);
Pd = input(2);
Hc = input(3);
Hd = input(4);

 % % SOC dynamics   

SoC_next = p.Ad*SoC_0 + (p.Bd*p.Delta_CT)*[Pc; Pd];
H_next = p.Ad_HSS*H_SS_0 + (p.Bd_HSS*p.Delta_CT)*[Hc; Hd];
end