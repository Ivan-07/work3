function [echo,Y_radar] = get_echo(angle,dis,S,v,c,F,M_T,M_R,rate,PRT,L)

tmp = 32.6 + 36.7*log10(dis);
alpha = 1/(dB_W(tmp)^(1/4));
% alpha = 1;
v_t = alpha*exp(1i*pi*sin(angle/180*pi)*(0:(M_T-1)));
v_r = alpha*exp(1i*pi*sin(angle/180*pi)*(0:(M_R-1)));
V = v_r.'*v_t/sqrt(rate);
tau = dis/(2*c);
f_d = v*F/(c*1e6);
% f_d = 0;
V = eye(M_T);
echo = zeros(M_R,PRT);
Y_radar = zeros(M_R,L);
for t = 1:L
    echo(:,t+tau) = exp(1i*2*pi*f_d*t)*V*S(:,t);
    Y_radar(:,t) = 1*V*S(:,t);
    %Y_radar(:,t)=  Y_radar(:,t)+ randn(size( Y_radar(:,t)));
end

end