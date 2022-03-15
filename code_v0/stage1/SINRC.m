function SINR_c = SINRC(L,h,w,p_r,g,S,sigma)

SINR_c = L*abs(h'*w)^2/(p_r*norm(g'*S)^2+L*sigma);