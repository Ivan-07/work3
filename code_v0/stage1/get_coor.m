function [coor_BS,coor_user] = get_coor(K,dis_Radar_BS,dis_BS_user)

coor_BS = zeros(K,2);
coor_user = zeros(K,2);
angle_diff = 180/K;
angle_res = (180-K*angle_diff)/2;
for k=1:K
    angle_Radar_BS = angle_res+(k-1)*angle_diff;
    coor_BS(k,1) = dis_Radar_BS*cos(angle_Radar_BS/180*pi);
    coor_BS(k,2) = dis_Radar_BS*sin(angle_Radar_BS/180*pi);
    
%     angle_user = round(rand*360);
    angle_BS_user = 90;
    coor_user(k,1) = coor_BS(k,1)+dis_BS_user*cos(angle_BS_user/180*pi);
    coor_user(k,2) = coor_BS(k,2)+dis_BS_user*sin(angle_BS_user/180*pi);
end



end