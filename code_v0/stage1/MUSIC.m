function angle_esti = MUSIC(Y_radar,L,M)
derad = pi/180; 
dd = 0.5;            % 阵元间距 
d=0:dd:(M-1)*dd;

X1 = Y_radar;
% 计算协方差矩阵
Rxx=X1*X1'/L;
% 特征值分解
[EV,D]=eig(Rxx);                   %特征值分解
EVA=diag(D)';                      %将特征值矩阵对角线提取并转为一行
[EVA,I]=sort(EVA);                 %将特征值排序 从小到大
EV=fliplr(EV(:,I));                % 对应特征矢量排序
                 
angle_esti = 0;
max_P = 0;
% 遍历每个角度，计算空间谱
for iang = 1:361
    angle(iang)=(iang-181)/2;
    phim=derad*angle(iang);
    a=exp(-1i*2*pi*d*sin(phim)).'; 
    En=EV(:,2:M);                   % 取矩阵的第M+1到N列组成噪声子空间
    Pmusic(iang)=1/(a'*En*En'*a);
    if abs(Pmusic(iang)) > max_P
        angle_esti = iang;
        max_P = abs(Pmusic(iang));
    end
end
Pmusic=abs(Pmusic);
Pmmax=max(Pmusic);
Pmusic=10*log10(Pmusic/Pmmax);            % 归一化处理
h=plot(angle,Pmusic);
set(h,'Linewidth',2);
xlabel('入射角/(degree)');
ylabel('空间谱/(dB)');
set(gca, 'XTick',[-90:30:90]);
grid on;

end