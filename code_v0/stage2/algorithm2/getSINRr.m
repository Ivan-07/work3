function SINR_r = getSINRr(u,a,S,w,G,sigma)

tmp1 = 0;tmp2=0;tmp3=0;

for l=1:L
    tmp1 = tmp1+abs(u(:,l)'*a(:,:,l)*S(:,l))^2;
    tmp2 = tmp2+abs(u(:,l)'*G*w)^2;
    tmp3 = tmp3+abs(u(:,l))^2;
end
SINR_r = tmp1/(tmp2+sigma*tmp3);

end