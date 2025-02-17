% Interpolate the state according to the following

% SoC(k+1) == Ad*SoC(k) + Bd*[Pc(k); Pd(k)]
% H(k+1) == Ad_HSS*H(k) + Bd_HSS*[Hc(k); Hd(k)]
%
function y = interpolate_pts(x_k,u_k,p)
y(:,1) = x_k;                     

    for t = 1:(p.Delta_DA/p.Delta_CT-1)
        y(:,t+1) = [p.Ad*y(1,t) + p.Bd*[u_k(1), u_k(2)]'*p.Delta_CT;
                p.Ad_HSS*y(2,t)+p.Bd_HSS*[u_k(3), u_k(4)]'*p.Delta_CT];
    end
end

% function y = interpolate_pts(x_k,u_k,delta,delta_low)
% y(:,1) = x_k;
% for t = 1:(delta/delta_low-1)
%     for i = 1:length(x_k)/2
%         y(2*i-1:2*i,t+1) = [x_k(2*i-1) + x_k(2*i)*t*delta_low + u_k(i)/2*(t*delta_low)^2;
%                 x_k(2*i)+u_k(i)*t*delta_low];
%     end
% %     y(:,t+1) = [x_k(1) + x_k(2)*t*delta_low + u_k(1)/2*(t*delta_low)^2;
% %                 x_k(2)+u_k(1)*t*delta_low;
% %                 x_k(3) + x_k(4)*t*delta_low + u_k(2)/2*(t*delta_low)^2;
% %                 x_k(4)+u_k(2)*t*delta_low;
% %                 x_k(5) + x_k(6)*t*delta_low + u_k(3)/2*(t*delta_low)^2;
% %                 x_k(6)+u_k(3)*t*delta_low;
% %                 ];
% end
% end
