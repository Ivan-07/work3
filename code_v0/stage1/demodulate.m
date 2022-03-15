function x_demodu = demodulate(Y_res,G,w,PRT)

x_demodu = zeros(1,PRT);
for t = 1:PRT
    cvx_begin quiet
    variable x complex
    minimize(sum_square_abs(Y_res(:,t)-G*w*x))
    cvx_end
    x_demodu(t) = x;
end


end