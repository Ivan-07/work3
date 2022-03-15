% function [tau_esti,f_d_esti] = estimate()
clc;clear all
K = 10;
M_T = 8;
M_R = 8;
L = 8;
M = 8;
N_T = 6;
sigma = dB_W(-100-30);
sigma_norm = 100;
PRT = 12;
v = 100;
c = 100;
F = 0.3e7;
max_iter = 6;

load('channel.mat','G_all')
G = G_all(:,:,1);
w = 0.25*(ones(N_T,1)+1i*ones(N_T,1));

tau_all = zeros(max_iter+1,1);
f_d_all = zeros(max_iter+1,1);
angle_real = 50;
dis_real = 800;
S = get_signal(M_T,L);
S = dftmtx(M);
[echo_radar,Y_radar] = get_echo(angle_real,dis_real,S,v,c,F,M_T,M_R,sigma/sigma_norm,PRT,L);
echo_comm = get_commSignal(PRT,G,w); 
echo = echo_radar; % + echo_comm;

angle_esti = MUSIC(Y_radar,L,M_R);
[tau_esti,f_d_esti] = estiTauAndFd(echo,S,PRT,L);
tau_all(1) = tau_esti;
f_d_all(1) = f_d_esti;
Y_radar = reconst_Yr(angle_esti,tau_esti,f_d_esti,c,M_T,M_R,sigma/sigma_norm,PRT,L,S);
Y_res = echo - Y_radar;
x_demodu = demodulate(Y_res,G,w,PRT);
Y_comm = G*w*x_demodu;
N = echo-Y_radar-Y_comm;

for i=1:max_iter
    Y_radar = Y_radar+N;
    [tau_esti,f_d_esti] = estiTauAndFd(Y_radar,S,PRT,L);
    tau_all(i+1) = tau_esti;
    f_d_all(i+1) = f_d_esti;
    Y_radar = reconst_Yr(angle_esti,tau_esti,f_d_esti,c,M_T,M_R,sigma/sigma_norm,PRT,L,S);
    Y_res = echo - Y_radar;
    x_demodu = demodulate(Y_res,G,w,PRT);
    Y_comm = G*w*x_demodu;
    N = echo-Y_radar-Y_comm;
end

[dis_real/(2*c) v*F/(c*1e6)]
[tau_all f_d_all]
