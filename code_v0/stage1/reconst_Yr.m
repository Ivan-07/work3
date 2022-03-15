function Y_R = reconst_Yr(angle,tau,f_d,c,M_T,M_R,rate,PRT,L,S)

dis = c*tau/2;
tmp = 32.6 + 36.7*log10(dis);
alpha = 1/(dB_W(tmp)^(1/4));
v_t = alpha*exp(1i*pi*sin(angle/180*pi)*(0:(M_T-1)));
v_r = alpha*exp(1i*pi*sin(angle/180*pi)*(0:(M_R-1)));
V = v_r*v_t'/sqrt(rate);

Y_R = zeros(M_R,PRT);
for t = 1:L
    Y_R(:,t+tau) = exp(1i*2*pi*f_d*t)*V*S(:,t);
end

end