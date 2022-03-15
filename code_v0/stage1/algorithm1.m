clear;clc;

K = 10;
M_T = 8;
M_R = 8;
L = 8;
M = 4;
N_T = 6;
Energy = zeros(K,1);
sigma = dB_W(-100-30);
sigma_norm = 10;
gamma_r = dB_W(20);
gamma_c = dB_W(20);

S_0 = get_signal(M_T,L);
% [a,G_all,h_all,g_all] = get_channel(M_R,L,M,M_T,S_0,N_T,K,sigma/sigma_norm);
% save('channel.mat','a','G_all','h_all','g_all')
load('channel.mat','a','G_all','h_all','g_all')

for k = 1:1
    G = G_all(:,:,k)*0.1;
    h = h_all(:,k)*0.1;
    g = g_all(:,k)*0.1;
    [u,w_k,p] = PDD(a,G,M,M_R,L,N_T,h,sigma_norm,gamma_r,gamma_c,g,S_0,M_T);
    fprintf("BS %d - E_c: %f \n sinr_c: %f \n ",k,norm(w_k)^2,SINRC(L,h,w_k,p^2,g,S_0,sigma_norm)/gamma_c)
    SINRR(M,L,u,a,G,w_k,p^2,sigma_norm)/gamma_r
    Energy(k) = L*norm(w_k,2)+p^2;
end