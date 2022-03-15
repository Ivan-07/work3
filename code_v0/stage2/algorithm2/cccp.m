function [S,w,v,u,x,y] = cccp(L,rho,a,G,x_pre,h,sigma,lambda1,lambda2,gamma_c,g,M_R,w,S,x,E_c,N_T,M_T,v,u,E_r)
eps = 1e-3;
iter_max = 30;
objValue_pre = 1;

for i=1:iter_max
    [v,x,y] = solve1(L,rho,u,a,G,x_pre,h,sigma,lambda1,lambda2,gamma_c,g,S,w);
    u = solve2(M_R,L,v,G,w,sigma,a,S);
    w = solve3(L,u,G,v,rho,h,x,lambda1,E_c,N_T);
    S = solve4(M_T,L,v,u,a,E_r,y,g,rho,lambda2);
    objValue = get_objValue(L,u,G,w,a,S,y,g,rho,lambda1,lambda2,v,sigma,x,h);
    fprintf("Inter %d - res: %f, %f\n",i, objValue, (objValue-objValue_pre));
    if(abs(objValue-objValue_pre) < eps)
        break;
    end
    objValue_pre = objValue;
end

end