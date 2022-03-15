function [u,w,p] = PDD(a,G,M,M_R,L,N_T,h,sigma,gamma_r,gamma_c,g,S_0,M_T)
[p,w,z,u,x,y] = init_X(M,M_R,L,N_T);
rho = 0.1*ones(3,1);
c = 0.5;
iter_max = 30;
eta = 1;
lambda1 = 0.1*ones(M,L);
lambda2 = 0.1*ones(M,L);
lambda3 = 0.1;

for i = 1:iter_max
    [p,w,z,u,x,y] = cccp(N_T,L,M,rho,u,a,x,y,G,z,h,sigma,lambda1,lambda2,lambda3,gamma_r,w,gamma_c,g,S_0,M_T);
%     SINRR(M,L,u,a,G,w,p^2,sigma)/gamma_r
    h_indicator = get_indicator(M,x,p,u,a,y,G,w,z,h,L);
    fprintf("iter(out): "+i+" "+h_indicator+"\n");
    if norm(h_indicator,'inf') <= eta
        for m=1:M
            u_m_L = u(:,:,m);
            a_m_L = a(:,:,m);
            for l=1:L
                lambda1(m,l) = lambda1(m,l)+1/rho(1)*(x(m,l)-p*u_m_L(:,l)'*a_m_L(:,l));
                lambda2(m,l) = lambda2(m,l)+1/rho(2)*(y(m,l)-u_m_L(:,l)'*G*w);
            end
        end
        lambda3 = lambda3+1/rho(3)*(z-h'*w);
    else
        rho = rho*c;
    end
    eta = 0.8*norm(h_indicator,'inf')
    if(norm(h,'inf') <= eps)
        break;
    end
end

end