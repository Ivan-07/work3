clear;clc;

K = 10;
M_T = 8;
M_R = 8;
L = 8;
M = 4;
N_T = 6;       
gamma_c = dB_W(20);
sigma = dB_W(-100-30);
sigma_norm = 10;
E_c = 1;
E_r = 1;

v = 50;
theta = 50;
dis = 200;

F = 1e8;
V = 3*1e8;
f_d = v*F/V;

% [~,G_all,h_all,g_all] = get_channel(M_R,L,M,M_T,S_0,N_T,K,sigma/sigma_norm);
% save('channel.mat','G_all','h_all','g_all')
load('channel.mat','G_all','h_all','g_all')

SINR_r = zeros(K,1);

for k = 1:K
    G = G_all(:,:,k);
    h = h_all(:,k);
    g = g_all(:,k);
    a = getA(theta,dis,f_d,M_T,M_R,L,sigma/sigma_norm);
    [u,w_k,S] = PDD(L,a,G,h,sigma,gamma_c,g,M_R,E_c,N_T,M_T,E_r);
    SINR_r(k) = getSINRr(u,a,S,w_k,G,sigma);
end