function ans = solve1_judge(gamma_c,y,x,sigma,L,x_pre)

ans = gamma_c*(norm(y)^2+L*sigma)+L*abs(x_pre)^2-2*L*real(conj(x_pre)*x);

end