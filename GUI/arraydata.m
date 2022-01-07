function [xs,xn,err] = arraydata(s,ndB,doa,coordinate_music,lambda)
%% 根据输入信号，干扰以及波达方向生成阵列数据快拍

Jamnum = length(doa(:,1));
elementnum = length(coordinate_music(:,1));
N = length(s);

% 加入误差矩阵
% 通道间增益不一致性小于1dB，均匀分布
% 通道间相位不一致性小于5°，均匀分布
gain_d = 1;
phase_d = 10;
delta_gain = gain_d*rand(elementnum,1)-gain_d/2;
delta_phase = phase_d*rand(elementnum,1)-phase_d/2;
err = 10.^(delta_gain/10).*exp(1i*pi*delta_phase./180);

%================================================================================
UintVector(1:Jamnum,:) = [sin(doa(:,1)).*cos(doa(:,2)),sin(doa(:,1)).*sin(doa(:,2)),cos(doa(:,1))];  % 干扰及信号方向向量
SteeringVector = exp(1i*2*pi/lambda*coordinate_music*UintVector');                                 % 干扰及信号导向矢量

xs = zeros(elementnum,N);
for i = 1:Jamnum
    xs = xs+SteeringVector(:,i)*s(i,:);% 形成干扰阵列数据
end
xn = wgn(elementnum,N,ndB);                                                % 各通道噪声，必须每个通道独立生成，不相关噪声
% xs = SteeringVector(:,Jamnum+1)*s;                                         % 形成信号阵列数据 

end