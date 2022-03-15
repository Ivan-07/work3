clear;clc;

K = 10;
M_T = 8;
M_R = 8;
L = 4;
M = 8;
N_T = 6;
Energy = zeros(K,1);
sigma = 1e-13;
gamma_r = 0.1;
gamma_c = 0.1;

S_0 = get_signal(M_T,L);
[a,G_all,h_all,g_all] = get_channel(M_R,L,M,M_T,S_0,N_T,K);

for k = 1:K
    G = G_all(:,:,k);
    h = h_all(:,k);
    g = g_all(:,k);
    [u,w_k,p] = PDD(a,G,M,M_R,L,N_T,h,sigma,gamma_r,gamma_c,g,S_0,M_T);
    Energy(k) = L*norm(w_k,2)+p^2;
end