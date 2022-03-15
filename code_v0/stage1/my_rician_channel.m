function [rician_channel] = my_rician_channel(m,n,dis,beta,sv1,sv2)
tmp = 32.6 + 36.7*log10(dis);
tmp = sqrt(dB_W(tmp));
L = 1/tmp;
g_random = 1/sqrt(2)*(randn(m ,n)+sqrt(-1)*(randn(m, n)));   
rician_channel = L*(sqrt(beta/(1+beta))*sv1*sv2'+sqrt(1/(1+beta))*g_random);
end