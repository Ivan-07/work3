function [a,G,h,g] = get_channel(M_R,L,M,M_T,S_0,N_T,K,rate)
angle_range = [225,315]+0;
dis = 1000;
tmp = 32.6 + 36.7*log10(dis);
alpha = 1/(dB_W(tmp)^(1/4));
a = zeros(M_R,L,M);
for m=1:M
    angle = angle_range(1) + (angle_range(2)-angle_range(1))/(M-1)*m;
    v_t = sqrt(1/(M_T))*alpha*exp(1i*pi*sin(angle/180*pi)*(0:(M_T-1)));
    v_r = sqrt(1/(M_R))*alpha*exp(1i*pi*sin(angle/180*pi)*(0:(M_R-1)));
    V = v_r*v_t'/sqrt(rate);
    for l=1:L
        a(:,l,m) = V*S_0(:,l);
    end
end

dis_Radar_BS = 100;
dis_BS_user = 20;

[coor_BS,coor_user] = get_coor(K,dis_Radar_BS,dis_BS_user);

angle_Radar_BS = zeros(1,K);
angle_BS_Radar = zeros(1,K);
angle_BS_user = zeros(1,K);
angle_Radar_user = zeros(1,K);
vec_Radar_BS = zeros(M_T,K);
vec_BS_Radar = zeros(N_T,K);
vec_BS_user = zeros(N_T,K);
vec_Radar_user = zeros(M_T,K);
dis_BS_user = zeros(1,K);
dis_Radar_user = zeros(1,K);
for k=1:K
    angle_Radar_BS(k) = mod(atan(coor_BS(2)/coor_BS(1))+pi,pi);
    vec_Radar_BS(:,k) = sqrt(1/(M_T))*exp(1i*pi*sin(angle_Radar_BS(k))*(0:(M_T-1)));
    angle_BS_Radar(k) = angle_Radar_BS(k)+pi;
    vec_BS_Radar(:,k) = sqrt(1/(N_T))*exp(1i*pi*sin(angle_BS_Radar(k))*(0:(N_T-1)));
    diff_BS_user = coor_user(k,:)-coor_BS(k,:);
    dis_BS_user(k) = norm(diff_BS_user);
    angle_BS_user(k) = atan(diff_BS_user(2)-diff_BS_user(1))+pi*(diff_BS_user(1)<0);
    vec_BS_user(:,k) = sqrt(1/(N_T))*exp(1i*pi*sin(angle_BS_user(k))*(0:(N_T-1)));
    dis_Radar_user(k) = norm(coor_user(k,:));
    angle_Radar_user(k) = mod(atan(coor_user(2)/coor_user(1))+pi,pi);
    vec_Radar_user(:,k) = sqrt(1/(M_T))*exp(1i*pi*sin(angle_Radar_user(k))*(0:(M_T-1)));
end


beta1 = dB_W(9);
beta2 = dB_W(3);

G = zeros(M_R,N_T,K);
h = zeros(N_T,K);
g = zeros(M_T,K);

for k=1:K
    G(:,:,k) = my_rician_channel(M_R,N_T,dis_Radar_BS,beta2,vec_Radar_BS(:,k),vec_BS_Radar(:,k))/sqrt(rate);
    h(:,k) = my_rician_channel(N_T,1,dis_BS_user(k),beta1,vec_BS_user(:,k),1)/sqrt(rate);
    g(:,k) = my_rician_channel(M_T,1,dis_Radar_user(k),beta2,vec_Radar_user(:,k),1)/sqrt(rate);
end



end