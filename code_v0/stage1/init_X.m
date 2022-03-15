function [p,w,z,u,x,y] = init_X(M,M_R,L,N_T)
% u = rand(M_R,L,M)+1i*rand(M_R,L,M);
% w = rand(N_T,1)+1i*rand(N_T,1);
% p = rand;
% z = rand+1i*rand;
% x = rand(M,L)+1i*rand(M,L);
% y = rand(M,L)+1i*rand(M,L);

u = 0.25*(ones(M_R,L,M)+1i*ones(M_R,L,M));
w = 0.25*(ones(N_T,1)+1i*ones(N_T,1));
p = 2;
z = 0.25+1i*0.25;
x = 0.25*(ones(M,L)+1i*ones(M,L));
y = 0.25*(ones(M,L)+1i*ones(M,L));
end