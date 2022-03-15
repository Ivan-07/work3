function objValue = get_objValue(L,w,p,rho,x,u,a,lambda1,y,G,z,h,lambda2,lambda3,M)

objValue = 0;
for m=1:M
    u_m = u(:,:,m);
    a_m = a(:,:,m);
    for l=1:L
        objValue = objValue+0.5/rho(1)*abs(x(m,l)-p*u_m(:,l)'*a_m(:,l)+rho(1)*lambda1(m,l))^2+...
            0.5/rho(2)*abs(y(m,l)-u(:,l)'*G*w+rho(2)*lambda2(m,l))^2;
    end
end
objValue = objValue+L*norm(w)^2+p^2+0.5/rho(3)*abs(z-h'*w+rho(3)*lambda3);
end