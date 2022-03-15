function h_indicator = get_indicator(M,x,p,u,a,y,G,w,z,h,L)

h_indicator = 0;

for m=1:M
    u_m_L = u(:,:,m);
    a_m_L = a(:,:,m);
    for l=1:L
        h_indicator = max(h_indicator,max(abs(x(m,l)-p*u_m_L(:,l)'*a_m_L(:,l)),abs(y(m,l)-u_m_L(:,l)'*G*w)));
    end
end

h_indicator = max(h_indicator,abs(z-h'*w));

end