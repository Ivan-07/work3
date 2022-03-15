function ans = solve1_judge(gamma_c,p,g,S,sigma,L,z_pre,z)

ans = gamma_c*(p^2*norm(g'*S)^2+L*sigma)+L*abs(z_pre)^2-2*L*real(conj(z_pre)*z);

end