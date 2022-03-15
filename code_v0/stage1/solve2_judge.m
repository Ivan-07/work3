function ans = solve2_judge(L,gamma_r,sigma,u_m,y_m,x_pre_m,x_m)

ans = 0;
for l=1:L
    ans = ans + gamma_r*(sigma*norm(u_m(:,l))^2+abs(y_m(l))^2)+...
        abs(x_pre_m(l))^2-2*real(conj(x_pre_m(l))*x_m(l));
end

end