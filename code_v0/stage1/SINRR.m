function SINR_r = SINRR(M,L,u,a,G,w,p_r,sigma)

SINR_r = zeros(M,1);
for m=1:M
    tmp1 = 0;
    tmp2 = 0;
    tmp3 = 0;
    u_m_L = u(:,:,m);
    a_m_L = a(:,:,m);
    for l=1:L
        tmp1 = tmp1+abs(u_m_L(:,l)'*a_m_L(:,l))^2;
        tmp2 = tmp2+abs(u_m_L(:,l)'*G*w)^2;
        tmp3 = tmp3+norm(u_m_L(:,l))^2;
    end
    SINR_r(m) = p_r*tmp1/(tmp2+sigma*tmp3);
end
    