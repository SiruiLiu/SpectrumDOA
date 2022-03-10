clear;
close all;
clc;

Fc = 1575420000;
C = 299792458;
Lambda = C/Fc;

IQData = load('IQData.mat');
Position = load('Position.mat');

IQData = IQData.SampledData;
Position = Position.Position;

% fp1 = fopen('.\sim_data\ch1_iq_data.txt', 'w+');
% fp2 = fopen('.\sim_data\ch2_iq_data.txt', 'w+');
% fp3 = fopen('.\sim_data\ch3_iq_data.txt', 'w+');
% fp4 = fopen('.\sim_data\ch4_iq_data.txt', 'w+');
% fp5 = fopen('.\sim_data\ch5_iq_data.txt', 'w+');
% fp6 = fopen('.\sim_data\ch6_iq_data.txt', 'w+');
% fp7 = fopen('.\sim_data\ch7_iq_data.txt', 'w+');
% fp8 = fopen('.\sim_data\ch8_iq_data.txt', 'w+');
% for i = 1:size(IQData, 2);
%     fprintf(fp1, '%f %f\r\n', real(IQData(1,i)), imag(IQData(1,i)));
%     fprintf(fp2, '%f %f\r\n', real(IQData(2,i)), imag(IQData(2,i)));
%     fprintf(fp3, '%f %f\r\n', real(IQData(3,i)), imag(IQData(3,i)));
%     fprintf(fp4, '%f %f\r\n', real(IQData(4,i)), imag(IQData(4,i)));
%     fprintf(fp5, '%f %f\r\n', real(IQData(5,i)), imag(IQData(5,i)));
%     fprintf(fp6, '%f %f\r\n', real(IQData(6,i)), imag(IQData(6,i)));
%     fprintf(fp7, '%f %f\r\n', real(IQData(7,i)), imag(IQData(7,i)));
%     fprintf(fp8, '%f %f\r\n', real(IQData(8,i)), imag(IQData(8,i)));
% end
% 
% fclose(fp1);
% fclose(fp2);
% fclose(fp3);
% fclose(fp4);
% fclose(fp5);
% fclose(fp6);
% fclose(fp7);
% fclose(fp8);

MUSIC_Simulation(IQData, Position, Lambda);