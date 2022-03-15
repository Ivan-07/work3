function [S,w,v,u,x,y] = init_X(M_R,L,N_T,M_T)
u = 0.5*(ones(M_R,L)+1i*ones(M_R,L));
w = 0.5*(ones(N_T,1)+1i*ones(N_T,1));
S = 0.25*(ones(M_T,L)+1i*ones(M_T,L));
v = 0.25+1i*0.25;
x = 0.25+1i*0.25;
y = 0.25*(ones(1,L)+1i*ones(1,L));
end