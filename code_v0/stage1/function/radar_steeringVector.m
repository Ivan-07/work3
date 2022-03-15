function a = radar_steeringVector(M,K)
angle_range = [135,255]+0;
a = zeros(M,K);
dis = 1000;
tmp = 32.6 + 36.7*log10(dis) - 0.7;
tmp = (dB_W(tmp))^(1/4);
L = 1/tmp;
for k=1:K
    angle = angle_range(1) + (angle_range(2)-angle_range(1))/(K-1)*k;
    a(:,k) = L*exp(1i*pi*sin(angle/180*pi)*(0:(M-1)));
end
end