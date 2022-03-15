function objValue = solve2_obj(L,a_m_L,rho,u_m,y_m,lambda1,x_m,lambda2,G,w,p,m)

objValue = 0;

for l=1:L
    objValue=objValue+abs(x_m(l)-p*u_m(:,l)'*a_m_L(:,l)+rho(1)*lambda1(m,l))^2+abs(y_m(l)-u_m(:,l)'*G*w+rho(2)*lambda2(m,l))^2;
end