function S_0 = get_signal(M_T,L)

S_0 = zeros(M_T,L);
for n=1:M_T
    for t = 1:L
        S_0(n,t) = sqrt(1/(M_T*L))*exp(1i*2*pi*n*(t-1)/L)*exp(1i*pi*(t-1)^2/L);
    end
end

end