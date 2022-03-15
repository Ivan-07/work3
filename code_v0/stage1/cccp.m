function [p,w,z,u,x,y] = cccp(N_T,L,M,rho,u,a,x,y,G,z,h,sigma,lambda1,lambda2,lambda3,gamma_r,w,gamma_c,g,S_0,M_T)
eps = 1e-3;
iter_max = 30;
objValue_pre = 1;

% objValues = [];

for i=1:iter_max
    [p,w,z] = solve1(N_T,L,M,rho,u,a,x,y,G,z,h,sigma,lambda1,lambda2,lambda3,gamma_c,g,S_0);
    [u,x,y] = solve2(M_T,L,M,gamma_r,sigma,G,w,x,a,rho,lambda2,p,lambda1,gamma_c);
    objValue = get_objValue(L,w,p,rho,x,u,a,lambda1,y,G,z,h,lambda2,lambda3,M);
    fprintf("Inter %d - res: %f, %f\n",i, objValue, (objValue-objValue_pre));
    if(abs(objValue-objValue_pre) < eps)
        break;
    end
    objValue_pre = objValue;
end

end