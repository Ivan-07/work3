function a = getA(theta,dis,f_d,M_T,M_R,L,rate)


tmp = 32.6 + 36.7*log10(dis) - 0.7;
alpha = 1/dB_W(tmp);
a = zeros(M_R,M_T,L);
v_t = sqrt(1/(M_T))*alpha*exp(1i*pi*sin(theta/180*pi)*(0:(M_T-1)));
v_r = sqrt(1/(M_R))*alpha*exp(1i*pi*sin(theta/180*pi)*(0:(M_R-1)));
V = v_r*v_t'/sqrt(rate);
for l=1:L
    a(:,:,l) = exp(1i*2*pi*f_d*l)*V;
end


end