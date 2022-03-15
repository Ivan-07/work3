function objValue = get_objValue(L,u,G,w,a,S,y,g,rho,lambda1,lambda2,v,sigma,x,h)

tmp1=0;tmp2=0;tmp3=0;
for l=1:L
    tmp1 = tmp1+abs(u(:,l)'*G*w)^2;
    tmp2 = tmp2+u(:,l)'*a(:,:,l)*S(:,l);
    tmp3 = tmp3+abs(y(l)-g'*S(:,l)+rho(2)*lambda2(l))^2;
end

objValue = abs(v)^2*(tmp1+sigma*norm(u)^2)-2*real(conj(v)*tmp2)+0.5/rho(1)*abs(x-h'*w+rho(1)*lambda1)^2+0.5/rho(2)*tmp3;

end