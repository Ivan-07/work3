function [tau_esti,f_d_esti] = estiTauAndFd(Y,S_0,PRT,L)

f_d_range = 1:500;
tau_range = 1:PRT-1-L;

f_d_esti = -1;
tau_esti = -1;

sum_conv = 0;
% for tau = tau_range
%     sum_tmp = 0;
%     for l=1:L
%         sum_tmp = sum_tmp + Y(:,tau+l)'*S_0(:,l);
%     end
%     if abs(sum_tmp) > sum_conv
%         tau_esti = tau;
%         sum_conv = abs(sum_tmp);
%     end
% end
% 
% sum_conv = 0;
% for f_d = f_d_range
%     sum_tmp = 0;
%     for l=1:L
%         sum_tmp = sum_tmp + Y(:,tau_esti+l)'*exp(-1i*2*pi*f_d*l)
%     end
%     if abs(sum_tmp) > sum_conv
%         f_d_esti = f_d;
%         sum_conv = abs(sum_tmp);
%     end
% end

sum_res = [];
for tau = tau_range
    for f_d = f_d_range
        sum_tmp = 0;
        for l=1:L
            sum_tmp = sum_tmp + sum(Y(:,tau+l).*S_0(:,l)', 'all')*exp(-1i*2*pi*f_d*(tau+l));
        end
        sum_res(tau, f_d) = sum_tmp;
        if abs(sum_tmp) > sum_conv
            tau_esti = tau;
            f_d_esti = f_d;
            sum_conv = abs(sum_tmp);
        end
    end
end

1
% [tau_esti f_d_esti]
