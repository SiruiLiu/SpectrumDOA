function [xs,xn,err] = arraydata(s,ndB,doa,coordinate_music,lambda)
%% ���������źţ������Լ����﷽�������������ݿ���

Jamnum = length(doa(:,1));
elementnum = length(coordinate_music(:,1));
N = length(s);

% ����������
% ͨ�������治һ����С��1dB�����ȷֲ�
% ͨ������λ��һ����С��5�㣬���ȷֲ�
gain_d = 1;
phase_d = 10;
delta_gain = gain_d*rand(elementnum,1)-gain_d/2;
delta_phase = phase_d*rand(elementnum,1)-phase_d/2;
err = 10.^(delta_gain/10).*exp(1i*pi*delta_phase./180);

%================================================================================
UintVector(1:Jamnum,:) = [sin(doa(:,1)).*cos(doa(:,2)),sin(doa(:,1)).*sin(doa(:,2)),cos(doa(:,1))];  % ���ż��źŷ�������
SteeringVector = exp(1i*2*pi/lambda*coordinate_music*UintVector');                                 % ���ż��źŵ���ʸ��

xs = zeros(elementnum,N);
for i = 1:Jamnum
    xs = xs+SteeringVector(:,i)*s(i,:);% �γɸ�����������
end
xn = wgn(elementnum,N,ndB);                                                % ��ͨ������������ÿ��ͨ���������ɣ����������
% xs = SteeringVector(:,Jamnum+1)*s;                                         % �γ��ź��������� 

end