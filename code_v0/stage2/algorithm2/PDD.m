function [u,w,S] = PDD(L,a,G,h,sigma,gamma_c,g,M_R,E_c,N_T,M_T,E_r)
[S,w,v,u,x,y] = init_X(M_R,L,N_T,M_T);
rho = 0.5*ones(2,1);
c = 0.8;
iter_max = 30;
eta = 1;
lambda1 = 0.1;
lambda2 = 0.1*ones(1,L);

for i = 1:iter_max
    [S,w,v,u,x,y] = cccp(L,rho,a,G,x,h,sigma,lambda1,lambda2,gamma_c,g,M_R,w,S,x,E_c,N_T,M_T,v,u,E_r);
    h_indicator = get_indicator(x,y,w,h,L,S,g);
    fprintf("iter(out): "+i+" "+h_indicator+"\n");
    if norm(h_indicator,'inf') <= eta
        for l=1:L
            lambda2(l) = lambda2(l)+1/rho(2)*(y(l)-g'*S(:,l));
        end
        lambda1 = lambda1+1/rho(1)*(x-h'*w);
    else
        rho = rho*c;
    end
    eta = 0.8*norm(h_indicator,'inf');
    if(norm(h,'inf') <= eps)
        break;
    end
end


end